from functools import wraps
import time


def timing_decorator(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        elapsed_time = time.time() - start_time
        print(f"Function {func.__name__} elapsed time: {elapsed_time:.4f} seconds")
        return result

    return wrapper
