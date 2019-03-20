# -*- coding: utf-8 -*-

def parse_vars(v):

    for key, value in v.items():
        if type(value) is int:
            print("integer : ", key, "=", value)
        if type(value) is float:
            print("float   : ", key, "=", value)
        if type(value) is list:
            print("list    : ", key, "=", value)
        if type(value) is bool:
            print("bool    : ", key, "=", value)

if __name__ == "__main__":

    a = 1
    b = 2.0
    c = [1, 2, 3]
    d = True

    parse_vars(vars())