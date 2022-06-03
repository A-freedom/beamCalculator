from beam import *
import unittest
from sympy import parse_expr
import json


class MyTestCase(unittest.TestCase):
    def test_beam(self):
        self.maxDiff = None
        with open("beam_test.json", "r") as read_file:
            testing_data = json.load(read_file)
            counter = 0
        for expect in testing_data['beams']:
            name = expect['name']
            # if name != 'Singer & Andre Problem 419':
            #     continue
            reactions = [Reaction(re[0], re[1]) for re in expect['reactions']]
            moments = [Moment(me[0], me[1]) for me in expect['moments']]
            loads = [Fun(parse_expr(lo[0]), lo[1]) for lo in expect['loads']]
            l = expect['l']
            actual = Beam(loads, reactions, moments, l, name=name).jsonDetails()
            with self.subTest(name):
                if not self.compareForces(expect['divided_loads'], actual['divided_loads']):
                    self.assertEqual(expect, actual)
                if not self.compareForces(expect['shear_stress'], actual['shear_stress']):
                    self.assertEqual(expect, actual)
                if not self.compareForces(expect['bending_moments'], actual['bending_moments']):
                    self.assertEqual(expect, actual)
                counter += 1
        print(f'{counter} beams have been successfully tested'.format(counter=counter))

    def compareForces(self, expect, actual):
        if len(expect) == len(actual):
            for i in range(len(expect)):
                if expect[i][1] != actual[i][1]: return False
                if parse_expr(expect[i][0]) != parse_expr(actual[i][0]): return False
        else:
            return False
        return True


if __name__ == '__main__':
    unittest.main()

# class MyTestCase(unittest.TestCase):
#     reactions = [Reaction(15, 1), Reaction(15, 3)]
#     loads = [Load(7.5 * x, [0, 2]), Load(15 - 7.5 * x, [2, 4])]
#     moments = []
#     beam = Beam(loads, reactions, moments, 4)
#
#     def test_domains(self):
#         self.assertEqual(self.beam.domains, [[0, 1], [1, 2], [2, 3], [3, 4]])
#
#     def test_divided_loads(self):
#         source_loads = [[[0, 1], -7.5 * x], [[1, 2], -7.5 * x - 7.5], [[2, 3], 7.5 * x - 15], [[3, 4], 7.5 * x - 7.5]]
#         self.assertEqualEquations(self.beam.divided_loads, source_loads)
#
#     def test_getShearEquations(self):
#         shear_source = [[[0, 1], -3.75 * x ** 2], [[1, 2], -3.75 * x ** 2 - 7.5 * x + 11.25],
#                         [[2, 3], 3.75 * x ** 2 - 15.0 * x + 0.e-124], [[3, 4], 3.75 * x ** 2 - 7.5 * x + 3.75]]
#         self.assertEqualEquations(self.beam.shear_stress, shear_source)
#
#     def test_getMomentEquations(self):
#         moment_source = [[[0, 1], -1.25 * x ** 3], [[1, 2], -1.25 * x ** 3 - 3.75 * x ** 2 + 11.25 * x - 1.25],
#                          [[2, 3], 1.25 * x ** 3 - 7.5 * x ** 2 + 1.89091402092252e-124 * x + 5.0],
#                          [[3, 4], 1.25 * x ** 3 - 3.75 * x ** 2 + 3.75 * x - 1.25]]
#
#         self.assertEqualEquations(self.beam.bending_moment, moment_source)
#
#     def assertEqualEquations(self, expect, actual):
#         source = [[i.interval, i.force] for i in expect]
#         if len(source) != len(actual):
#             self.fail('expect {len1} equations but actual should be {len2} equations'.format(len1=len(source),
#                                                                                              len2=len(actual)))
#
#         for i in range(len(source)):
#             # the +1-1 is used to force the sympy to re calculate the equation and got out of terms like 0.e-124
#             a = str(source[i][1] + 1 - 1)
#             b = str(actual[i][1] + 1 - 1)
#             if source[i][0] != actual[i][0] or a != b:
#                 self.fail('Item of index {i} is different \n {h} != {k}'.format(i=i, h=source[i], k=actual[i]))
#
#
# if __name__ == '__main__':
#     unittest.main()
