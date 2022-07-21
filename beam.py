import numpy as np
import plotly.graph_objects as go

from plotly.subplots import make_subplots
from functools import reduce
import json
from sympy import Eq, integrate, symbols, Abs, expand, solve, lambdify, diff
from sympy.abc import x


class Beam:
    simple = 100

    def __init__(self, loads, reactions, moments, l, name='', simple=100, supports=None, I=1, E=0):
        if supports is None:
            supports = []
        self.I = I
        self.E = E
        self.supports = supports
        self.simple = simple
        self.name = name
        self.domains = []
        self.divided_loads = []
        self.shear_stress = []
        self.bending_moment = []
        self.deflectionSlope = []
        self.deflections = []
        self.C1 = 0
        self.C2 = 0
        self.loads = loads
        self.reactions = reactions
        self.moments = moments
        self.l = l
        self.calculateReactions()
        self.setDomains()
        self.setDivided_loads()
        self.calculateShearForce()
        self.calculateBendingMoment()
        self.calculateDeflection()

    def setDomains(self):
        intervals_list = [0, self.l]
        for i in self.loads:
            intervals_list += i.interval
        for i in self.reactions:
            intervals_list.append(i.offset)
        for i in self.moments:
            intervals_list.append(i.offset)
        self.domains = []
        intervals_list.sort()
        intervals_list = list(dict.fromkeys(intervals_list))
        for i in range(len(intervals_list) - 1):
            self.domains.append([intervals_list[i], intervals_list[i + 1]])

    def setDivided_loads(self):
        for domain in self.domains:
            tem_load = [Abs(0)]
            for lo in self.loads:
                a1 = lo.interval[0]
                a2 = lo.interval[1]
                b1 = domain[0]
                b2 = domain[1]
                if not (a1 >= b2 or b1 >= a2):
                    tem_load.append(- lo.force.subs(x, x + b1 - a1))
            self.divided_loads.append(RawLoad(reduce(lambda a, b: a + b, tem_load), domain))

    def calculateReactions(self):
        equations = []
        anonP = ()
        anonM = ()
        Cm = 1
        Cp = 1
        for i in self.supports:
            symbol = symbols('P{c}'.format(c=Cp))
            self.reactions.append(PointedLoad(symbol, i.offset))
            anonP += (symbol,)
            Cp += 1
            if i.type == 'fix':
                symbol = symbols('M{c}'.format(c=Cm))
                self.moments.append(Moment(symbol, i.offset))
                anonM = (symbol,)
                Cm += 1

        for g in self.supports:
            summation = Abs(0)
            for i in self.moments:
                summation += i.force
            for i in self.reactions:
                summation += i.force * (g.offset - i.offset)
            for i in self.loads:
                a = i.interval[0]
                b = i.interval[1]
                A = integrate(-i.force, (x, a, b))
                Ax = integrate(-i.force * x, (x, a, b))
                M = A * g.offset - Ax
                summation += M
            equations.append(Eq(summation, 0))
        summationOfy = Abs(0)
        for i in self.loads:
            summationOfy += integrate(-i.force, (x, i.interval[0], i.interval[1]))
        for i in self.reactions:
            summationOfy += i.force
        equations.append(Eq(summationOfy, 0))
        _solve = solve(equations, anonP + anonM)
        print(_solve)
        self.reactions = [PointedLoad(i.force.subs(_solve), i.offset) for i in self.reactions]
        self.moments = [Moment(i.force.subs(_solve), i.offset) for i in self.moments]

    def reactionsAt(self, offset):
        re_ls = [0]
        for i in self.reactions:
            if i.offset == offset:
                re_ls.append(i.force)
        return reduce(lambda a, b: a + b, re_ls)

    def calculateShearForce(self):
        shear_stress = []
        for i in self.divided_loads:
            c = self.reactionsAt(i.interval[0])
            if len(shear_stress) != 0:
                off = shear_stress[-1].interval
                c += shear_stress[-1].force.evalf(subs={x: off[1] - off[0]})
            shear_stress.append(RawLoad(expand(integrate(i.force, x) + c), i.interval))
        self.shear_stress = shear_stress

    def momentAt(self, offset):
        re_ls = [0]
        for i in self.moments:
            if i.offset == offset:
                re_ls.append(i.force)
        return reduce(lambda a, b: a + b, re_ls)

    def calculateBendingMoment(self):
        bending_moment = []
        for i in self.shear_stress:
            c = self.momentAt(i.interval[0])
            if len(bending_moment) != 0:
                off = bending_moment[-1].interval
                c += bending_moment[-1].force.evalf(subs={x: off[1] - off[0]})
            bending_moment.append(RawLoad(expand(integrate(i.force, x) + c), i.interval))
        self.bending_moment = bending_moment

    def getFunctions(self, _list):
        return [lambdify(x, exp.force) for exp in _list]

    def printDetails(self):
        def export(eq):
            for i in eq:
                ex = self.findMixAndMin(i.force, i.interval)
                print(i.interval, '-->', i.force, )
                print('maximum = ', ex[0].y, ', at x = ', ex[0].x)
                print('minimum = ', ex[1].y, ', at x = ', ex[1].x)

        # print equations
        print("reactions ::")
        [print(i.offset, '-->', i.force) for i in self.reactions]
        print("moments ::")
        [print(i.offset, '-->', i.force) for i in self.moments]
        print("divided loads ::")
        [print(i.interval, '-->', i.force) for i in self.divided_loads]
        print("shear stress equations ::")
        export(self.shear_stress)
        print("bending moment equations ::")
        export(self.bending_moment)
        print("deflection slope ::")
        export(self.deflectionSlope)
        print("deflection ::")
        export(self.deflections)

    def jsonDetails(self):
        data = {"name": self.name,
                "l": self.l,
                "loads": [[str(i.force), i.interval] for i in self.loads],
                "reactions": [[i.force, i.offset] for i in self.reactions],
                "moments": [[i.force, i.offset] for i in self.moments],
                "divided_loads": [[str(i.force), i.interval] for i in self.divided_loads],
                "shear_stress": [[str(i.force), i.interval] for i in self.shear_stress],
                'bending_moments': [[str(i.force), i.interval] for i in self.bending_moment]
                }
        return data

    def displayPlots(self):
        def plot(row, col, color, name):
            fig.add_trace(
                go.Scatter(
                    x=xx,
                    y=yy,
                    fill='tonexty',
                    mode='lines',
                    name=name,
                    line_color=color,
                    line=dict(width=3)),
                row=row, col=col
            )
            fig.add_trace(
                go.Scatter(
                    x=[0, self.l],
                    y=[0, 0],
                    fill='tonexty',
                    mode='lines',
                    line_color=color,
                    line=dict(width=3)),
                row=row, col=col
            )

        fig = make_subplots(
            subplot_titles=("SHEAR STRESS", "BENDING MOMENTS", "THE SLOP", "DEFLECTIONS"),
            rows=2, cols=2,
            # shared_xaxes=True,
            print_grid=True,
            # column_widths=[200],
            # row_heights=[200]*4,
            vertical_spacing=0.03,
            # specs=[[{"type": "scatter"},{"type": "scatter"}],[{"type": "scatter"},{"type": "scatter"}]]
        )
        # shear
        yy = [0]
        shear_stress_fast = self.getFunctions(self.shear_stress)
        for i in range(len(self.domains)):
            for d in np.linspace(0, int(self.domains[i][1] - self.domains[i][0]),
                                 int((self.domains[i][1] - self.domains[i][0]) * self.simple)):
                yy.append(shear_stress_fast[i](d))
        yy.append(yy[-1] + self.reactionsAt(self.l).__float__())
        xx = np.linspace(0, self.l, len(yy))
        plot(1, 1, '#4D4DFF', 'shear stress')

        # bending
        yy = [0]
        bending_moment_fast = self.getFunctions(self.bending_moment)
        for i in range(len(self.domains)):
            for d in np.linspace(0, int(self.domains[i][1] - self.domains[i][0]),
                                 int((self.domains[i][1] - self.domains[i][0]) * self.simple)):
                yy.append(bending_moment_fast[i](d))
        yy.append(yy[-1] + self.momentAt(self.l).__float__())
        xx = np.linspace(0, self.l, len(yy))
        plot(1, 2, '#D22730', 'bending moments')


        # deflection slope
        yy = []
        slop_fast = self.getFunctions(self.deflectionSlope)
        for i in range(len(self.domains)):
            for d in np.linspace(0, int(self.domains[i][1] - self.domains[i][0]),
                                 int((self.domains[i][1] - self.domains[i][0]) * self.simple)):
                yy.append(slop_fast[i](d))
        xx = np.linspace(0, self.l, len(yy))
        plot(2, 1, '#FFAD00', 'deflection slop')

        # deflection
        yy = []
        deflections_fast = self.getFunctions(self.deflections)
        for i in range(len(self.domains)):
            for d in np.linspace(0, int(self.domains[i][1] - self.domains[i][0]),
                                 int((self.domains[i][1] - self.domains[i][0]) * self.simple)):
                yy.append(deflections_fast[i](d))
        xx = np.linspace(0, self.l, len(yy))
        plot(2, 2, 'green', 'deflection')

        fig.update_layout(
            showlegend=False,
            title_text="BEAM TEST",
            template='plotly_white'
        )
        fig.show()

    def save(self):
        def write_json(new_data, filename='beam_test.json'):
            with open(filename, 'r+') as file:
                # First we load existing data into a dict.
                file_data = json.load(file)
                # Join new_data with file_data inside emp_details
                file_data["beams"].append(new_data)
                # Sets file's current position at offset.
                file.seek(0)
                # convert back to json.
                json.dump(file_data, file)

        write_json(self.jsonDetails())

    def calculateDeflection(self):
        cost = symbols('C1:{c}'.format(c=1 + len(self.bending_moment) * 2))
        counter = 0
        for i in self.bending_moment:
            self.deflectionSlope.append(RawLoad(integrate(i.force, x) + cost[counter], i.interval))
            counter += 1

        for i in self.deflectionSlope:
            self.deflections.append(RawLoad(integrate(i.force, x) + cost[counter], i.interval))
            counter += 1

        equation = []
        for i in self.supports:
            if i.type == 'pin':
                for d in self.deflections:
                    if d.interval[0] == i.offset:
                        equation.append(Eq(0, d.force.replace(x, 0)))
                if self.deflections[-1].interval[1] == i.offset:
                    equation.append(Eq(0, self.deflections[-1].force.replace(x, self.deflections[-1].interval[1] -
                                                                             self.deflections[-1].interval[0])))
            elif i.type == 'fix':
                for d in self.deflections:
                    if d.interval[0] == i.offset:
                        equation.append(Eq(0, d.force.replace(x, 0)))
                    elif d.interval[1] == i.offset:
                        equation.append(Eq(0, d.force.replace(x, d.interval[1] - d.interval[0])))
                for d in self.deflectionSlope:
                    if d.interval[0] == i.offset:
                        equation.append(Eq(0, d.force.replace(x, 0)))
                    elif d.interval[1] == i.offset:
                        equation.append(Eq(0, d.force.replace(x, d.interval[1] - d.interval[0])))
        for i in range(1, len(self.deflectionSlope)):
            d = self.deflectionSlope[i - 1].interval[1] - self.deflectionSlope[i - 1].interval[0]
            equation.append(Eq(self.deflectionSlope[i - 1].force.replace(x, d),
                               self.deflectionSlope[i].force.replace(x, 0)))
            equation.append(Eq(self.deflections[i - 1].force.replace(x, d),
                               self.deflections[i].force.replace(x, 0)))

        _solve = solve(equation, cost)
        print(_solve)
        self.deflectionSlope = [
            RawLoad(i.force.subs(_solve) / (self.I * self.E), i.interval) for i in
            self.deflectionSlope]
        self.deflections = [RawLoad(i.force.subs(_solve) / (self.I * self.E), i.interval)
                            for i in self.deflections]

    def findMixAndMin(self, function, interval):
        firstD = diff(function, x)
        ext = solve(firstD, x)
        points = [Point(i - interval[0], function.subs(x, i - interval[0])) for i in interval]
        for i in ext:
            d = i
            if not i.is_real:
                d = i.args[0]
            if 0 < d < interval[1] - interval[0]:
                points.append(Point(d, function.subs(x, d)))

        _max = points[0]
        _min = points[0]
        for i in points:
            if i.y > _max.y:
                _max = i
            if i.y < _min.y:
                _min = i
        return [_max, _min]


class RawLoad:
    def __init__(self, load, interval):
        self.force = (Abs(0) + load)
        self.interval = interval


class Moment:
    def __init__(self, force, offset):
        self.force = force + Abs(0)
        self.offset = offset


class PointedLoad:
    def __init__(self, force, offset):
        self.force = force + Abs(0)
        self.offset = offset


class Support:
    def __init__(self, _type, offset):
        self.type = _type
        self.offset = offset + Abs(0)


class DistributedLoad:
    def __init__(self, start_offset, end_offset, start_value, end_value):
        self.start_offset = start_offset
        self.end_offset = end_offset
        self.start_value = start_value
        self.end_value = end_value
        self.force = (Abs(0) + (start_value + x * ((end_value - start_value) / (end_offset - start_offset))))
        self.interval = [start_offset, end_offset]


class Point:
    def __init__(self, _x, _y):
        self.x = _x
        self.y = _y
