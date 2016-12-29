from collections import defaultdict

def f(ab="123", bc="456", *args):
    print type(args)
    print len(args)
    for arg in args:
        print arg

f("bbb", "aaa", ["1", "2", "3"])