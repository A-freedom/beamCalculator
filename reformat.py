import re
from decimal import Decimal
from sympy import Float
from sympy.abc import x

# future me sorry for all the suffering you facing while reading this code try to learn from my mistake & write
# better code ***/
def reformat(eq):
    st = str(eq).split(' ')
    for i in range(len(st)):
        r = re.findall(r'-?[\d.]+(?:e-?\d+)?', st[i])
        for d in r:
            if Float(d) == 0:
                st[i] = st[i].replace(d, '0')
                continue
            if d == '+' or d == '-':
                continue
            if not (-0.1 < Float(d) < 0.1):
                if int(Float(d)) == Float(d):
                    st[i] = st[i].replace(d, str(int(Float(d))))
                    continue
                st[i] = st[i].replace(d, str(Float(d).round(2)))
                continue
            st[i] = st[i].replace(d, '%.2E' % Decimal(d))

    return ' '.join(st)
