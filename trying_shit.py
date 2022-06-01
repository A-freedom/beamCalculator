from beam import *
from sympy.parsing.sympy_parser import parse_expr

x = smp.symbols('x', real=True)
import json

if __name__ == '__main__':
    with open("beam_test.json", "r") as read_file:
        x = json.load(read_file)
    for i in x.get('beams'):
        reactions = [Reaction(re[0], re[1]) for re in i.get('reactions')]
        moments = [Moment(me[0], me[1]) for me in i.get('moments')]
        loads = [Load(parse_expr(lo[0]), lo[1]) for lo in i.get('loads')]
        l = i.get('l')
        beam = Beam(loads, reactions, moments, l)
        print(beam.jsonDetails())
        act = json.loads(beam.jsonDetails())
