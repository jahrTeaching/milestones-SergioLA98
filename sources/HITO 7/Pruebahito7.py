from psutil import cpu_count
from pyopencl.tools import get_test_platform_and_devices



def prueba():
    print(f'thread count per core: {cpu_count() // cpu_count(logical=False)}')


#prueba()
get_test_platform_and_devices()
