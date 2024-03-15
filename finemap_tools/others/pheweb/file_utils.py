from .utils import PheWebError, get_phenolist, chrom_order

# from . import conf
from . import parse_utils

import io
import os
import csv
from contextlib import contextmanager

# import json
import gzip

# import datetime
from boltons.fileutils import AtomicSaver, mkdir_p

# import pysam
import itertools, random
from pathlib import Path
from typing import List, Callable, Dict, Union, Iterator, Optional, Any

csv.register_dialect(
    "pheweb-internal-dialect",
    delimiter="\t",
    doublequote=False,
    escapechar="\\",
    lineterminator="\n",
    quotechar='"',
    skipinitialspace=False,
    strict=True,
)


def get_filepath(kind: str, *, must_exist: bool = True) -> str:
    if kind not in _single_filepaths:
        raise Exception("Unknown kind of filepath: {}".format(repr(kind)))
    filepath: str = _single_filepaths[kind]()
    if must_exist and not os.path.exists(filepath):
        raise PheWebError(
            "Filepath {} of kind {} was requested but doesn't exist".format(
                filepath, kind
            )
        )
    return filepath


@contextmanager
def VariantFileWriter(
    filepath: str, allow_extra_fields: bool = False, use_gzip: bool = True
):
    """
    Writes variants (represented by dictionaries) to an internal file.

        with VariantFileWriter('a.tsv') as writer:
            writer.write({'chrom': '2', 'pos': 47, ...})

    Each variant/association/hit/loci written must have a subset of the keys of the first one.
    """
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    if use_gzip:
        with AtomicSaver(
            filepath,
            text_mode=False,
            # part_file=part_file,
            overwrite_part=True,
            rm_part_on_exc=False,
        ) as f:
            with gzip.open(f, "wt", compresslevel=2) as f_gzip:
                yield _vfw(f_gzip, allow_extra_fields, filepath)
    else:
        with AtomicSaver(
            filepath,
            text_mode=True,
            # part_file=part_file,
            overwrite_part=True,
            rm_part_on_exc=False,
        ) as f:
            yield _vfw(f, allow_extra_fields, filepath)


class _vfw:
    def __init__(self, f, allow_extra_fields: bool, filepath: str):
        self._f = f
        self._allow_extra_fields = allow_extra_fields
        self._filepath = filepath

    def write(self, variant: Dict[str, Any]) -> None:
        if not hasattr(self, "_writer"):
            fields: List[str] = []
            for field in parse_utils.fields:
                if field in variant:
                    fields.append(field)
            extra_fields = list(set(variant.keys()) - set(fields))
            if extra_fields:
                if not self._allow_extra_fields:
                    raise PheWebError(
                        "ERROR: found unexpected fields {!r} among the expected fields {!r} while writing {!r}.".format(
                            extra_fields, fields, self._filepath
                        )
                    )
                fields += extra_fields
            self._writer = csv.DictWriter(
                self._f, fieldnames=fields, dialect="pheweb-internal-dialect"
            )
            self._writer.writeheader()
        self._writer.writerow(variant)

    def write_all(self, variants: Iterator[Dict[str, Any]]) -> None:
        for v in variants:
            self.write(v)


def make_basedir(path: Union[str, Path]) -> None:
    mkdir_p(os.path.dirname(path))


## reader


@contextmanager
def VariantFileReader(
    filepath: Union[str, Path], only_per_variant_fields: bool = False
):
    """
    Reads variants (as dictionaries) from an internal file.  Iterable.  Exposes `.fields`.

        with VariantFileReader('a.tsv') as reader:
            print(reader.fields)
            for variant in reader:
                print(variant)
    """
    with read_maybe_gzip(filepath) as f:
        reader: Iterator[List[str]] = csv.reader(f, dialect="pheweb-internal-dialect")
        try:
            fields = next(reader)
        except StopIteration:
            raise PheWebError("It looks like the file {} is empty".format(filepath))
        if fields[0].startswith(
            "#"
        ):  # This won't happen in normal use but it's convenient for temporary internal re-routing
            fields[0] = fields[0][1:]
        # for field in fields:
        #     assert (
        #         field in parse_utils.per_variant_fields
        #         or field in parse_utils.per_assoc_fields
        #     ), field

        if only_per_variant_fields:
            yield _vfr_only_per_variant_fields(fields, reader)
        else:
            yield _vfr(fields, reader)


@contextmanager
def read_maybe_gzip(filepath: Union[str, Path]):
    if isinstance(filepath, Path):
        filepath = str(filepath)
    is_gzip = False
    with open(filepath, "rb", buffering=0) as raw_f:  # no need for buffers
        if raw_f.read(3) == b"\x1f\x8b\x08":
            is_gzip = True
    if is_gzip:
        with read_gzip(filepath) as f:
            yield f
    else:
        with open(filepath, "rt", buffering=2**18) as f:  # 256KB buffer
            yield f


@contextmanager
def read_gzip(filepath):  # mypy doesn't like it
    # hopefully faster than `gzip.open(filepath, 'rt')` -- TODO: find out whether it is
    with gzip.GzipFile(
        filepath, "rb"
    ) as f:  # leave in binary mode (default), let TextIOWrapper decode
        with io.BufferedReader(f, buffer_size=2**18) as g:  # 256KB buffer
            with io.TextIOWrapper(g) as h:  # bytes -> unicode
                yield h


class _vfr_only_per_variant_fields:
    def __init__(self, fields: List[str], reader: Iterator[List[str]]):
        self._all_fields = fields
        self._extractors = [
            (parse_utils.reader_for_field[field], field, colidx)
            for colidx, field in enumerate(fields)
            if field in parse_utils.per_variant_fields
        ]
        self.fields = [e[1] for e in self._extractors]
        self._reader = reader

    def __iter__(self) -> Iterator[Dict[str, Any]]:
        return self._get_variants()

    def _get_variants(self) -> Iterator[Dict[str, Any]]:
        for unparsed_variant in self._reader:
            assert len(unparsed_variant) == len(self._all_fields), (
                unparsed_variant,
                self._all_fields,
            )
            variant = {
                field: parser(unparsed_variant[colidx])
                for parser, field, colidx in self._extractors
            }
            yield variant


class _vfr:
    def __init__(self, fields: List[str], reader: Iterator[List[str]]):
        self.fields = fields
        self._reader = reader

    def __iter__(self) -> Iterator[Dict[str, Any]]:
        return self._get_variants()

    def _get_variants(self) -> Iterator[Dict[str, Any]]:

        # print(self.fields)
        parsers: List[Callable[[str], Any]] = [
            parse_utils.reader_for_field.get(field, None) for field in self.fields
        ]
        for unparsed_variant in self._reader:
            assert len(unparsed_variant) == len(self.fields), (
                unparsed_variant,
                self.fields,
            )

            variant = {
                field: parser(value) if parser else value
                for parser, field, value in zip(parsers, self.fields, unparsed_variant)
            }
            # variant = {
            #     field: value for field, value in zip(self.fields, unparsed_variant)
            # }

            yield variant
