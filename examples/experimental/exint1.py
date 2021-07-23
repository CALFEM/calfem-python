# -*- coding: utf-8 -*-

import calfem.intvis as iv

# --- Parametrar som kan ändras

a = 1
b_slider = 2.0  
c_list = [1, 2, 3]
d_check = True
f_param = 42.0
g_float = 84.0
g_int = 34

# --- Redigera parametrar interaktivt

iv.edit_params(vars())

def test():
    c = 1
    d = 2
    iv.edit_params(vars())

test()

# --- Redigera geometri

#iv.edit_geometry(g)

