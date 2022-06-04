import numpy as np
import sympy as smp
import matplotlib.pyplot as plt
from functools import reduce
import json

from sympy import Eq
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
        self.loads = loads
        self.reactions = reactions
        self.moments = moments
        self.l = l
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
            tem_load = [smp.Abs(0)]
            for lo in self.loads:
                a1 = lo.interval[0]
                a2 = lo.interval[1]
                b1 = domain[0]
                b2 = domain[1]
                if not (a1 >= b2 or b1 >= a2):
                    tem_load.append(- lo.force.subs(x, x + b1 - a1))
            self.divided_loads.append(Fun(reduce(lambda a, b: a + b, tem_load), domain))

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
            shear_stress.append(Fun(smp.expand(smp.integrate(i.force, x) + c), i.interval))
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
            bending_moment.append(Fun(smp.expand(smp.integrate(i.force, x) + c), i.interval))
        self.bending_moment = bending_moment

    def getFunctions(self, _list):
        return [smp.lambdify(x, exp.force) for exp in _list]

    def printDetails(self):
        # print equations
        print("divided loads ::")
        [print(i.interval, '-->', i.force) for i in self.divided_loads]
        print("shear stress equations ::")
        [print(i.interval, '-->', i.force) for i in self.shear_stress]
        print("bending moment equations ::")
        [print(i.interval, '-->', i.force) for i in self.bending_moment]
        print("deflection slope ::")
        [print(i.interval, '-->', i.force) for i in self.deflectionSlope]
        print("deflection ::")
        [print(i.interval, '-->', i.force) for i in self.deflections]

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
        plt.style.use('fast')
        fig, ax = plt.subplots(1, 4, num=self.name)
        # shear
        yy = [0]
        shear_stress_fast = self.getFunctions(self.shear_stress)
        for i in range(len(self.domains)):
            for d in np.linspace(0, self.domains[i][1] - self.domains[i][0],
                                 (self.domains[i][1] - self.domains[i][0]) * self.simple):
                yy.append(shear_stress_fast[i](d))
        yy.append(yy[-1] + self.reactionsAt(self.l))
        xx = np.linspace(0, self.l, len(yy))
        ax[0].stackplot(xx, yy, color='blue', alpha=0.3)
        ax[0].plot(xx, yy, color='blue', linewidth=2.5)
        ax[0].grid(visible=True, which='major', axis='both')
        ax[0].axhline(y=0, color='k')
        ax[0].set_title("shear stress")

        # bending
        yy = [0]
        bending_moment_fast = self.getFunctions(self.bending_moment)
        for i in range(len(self.domains)):
            for d in np.linspace(0, self.domains[i][1] - self.domains[i][0],
                                 (self.domains[i][1] - self.domains[i][0]) * self.simple):
                yy.append(bending_moment_fast[i](d))
        yy.append(yy[-1] + self.momentAt(self.l))
        xx = np.linspace(0, self.l, len(yy))
        ax[1].stackplot(xx, yy, color='green', alpha=0.3)
        ax[1].plot(xx, yy, color='green', linewidth=2.5)
        ax[1].grid(visible=True, which='major', axis='both')
        ax[1].axhline(y=0, color='k')
        ax[1].set_title("bending moment")

        # deflection slope
        yy = []
        deflectionSlope_fast = self.getFunctions(self.deflectionSlope)
        for i in range(len(self.domains)):
            for d in np.linspace(0, self.domains[i][1] - self.domains[i][0],
                                 (self.domains[i][1] - self.domains[i][0]) * self.simple):
                yy.append(deflectionSlope_fast[i](d))
        xx = np.linspace(0, self.l, len(yy))
        ax[2].stackplot(xx, yy, color='green', alpha=0.3)
        ax[2].plot(xx, yy, color='green', linewidth=2.5)
        ax[2].grid(visible=True, which='major', axis='both')
        ax[2].axhline(y=0, color='k')
        ax[2].set_title("deflection slope")

        # deflection
        yy = []
        deflection_fast = self.getFunctions(self.deflections)
        for i in range(len(self.domains)):
            for d in np.linspace(0, self.domains[i][1] - self.domains[i][0],
                                 (self.domains[i][1] - self.domains[i][0]) * self.simple):
                yy.append(deflection_fast[i](d))
        xx = np.linspace(0, self.l, len(yy))
        ax[3].stackplot(xx, yy, color='green', alpha=0.3)
        ax[3].plot(xx, yy, color='green', linewidth=2.5)
        ax[3].grid(visible=True, which='major', axis='both')
        ax[3].axhline(y=0, color='k')
        ax[3].set_title("deflection")
        # show all
        plt.show()

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
        cost = smp.symbols('C1:{c}'.format(c=1 + len(self.bending_moment) * 2))
        counter = 0
        for i in self.bending_moment:
            self.deflectionSlope.append(Fun(smp.integrate(i.force, x) + cost[counter], i.interval))
            counter += 1
        for i in self.deflectionSlope:
            self.deflections.append(Fun(smp.integrate(i.force, x) + cost[counter], i.interval))
            counter += 1

        equation = []
        d_counter = 1
        for support in self.supports:
            if support.type == 'pin':
                for deflection in self.deflections:
                    if deflection.interval[0] == support.offset:
                        eq = Eq(deflection.force.replace(x, 0), 0)
                        print('d{c}(0) = 0 --> {eq} = 0'.format(c=d_counter, eq=eq.args[0]))
                        d_counter += 1
                        equation.append(eq)
                    elif deflection.interval[1] == support.offset:
                        eq = Eq(deflection.force.replace(x, deflection.interval[1] - deflection.interval[0]), 0)
                        print('d{c}({d}) = 0 --> {eq1} = {eq2}'.format(c=d_counter, eq1=eq.args[0], eq2=eq.args[1],
                                                                       d=deflection.interval[1] - deflection.interval[
                                                                           0]))
                        d_counter += 1
                        equation.append(eq)
        s_conter = 0
        d_counter = 1
        for i in range(1, len(self.deflections)):
            nextDeflection = self.deflections[i]
            deflection = self.deflections[i - 1]
            nextSlop = self.deflections[i]
            slop = self.deflections[i - 1]
            eq = Eq(deflection.force.replace(x, deflection.interval[1] - deflection.interval[0]),
                    nextDeflection.force.replace(x, 0))
            print('d{c}({d}) = d{cp}(0) --> {eq1} = {eq2}'.format(c=d_counter, eq1=eq.args[0], eq2=eq.args[1],
                                                                  cp=d_counter + 1,
                                                                  d=deflection.interval[1] - deflection.interval[
                                                                      0]))
            equation.append(eq)
            eq = Eq(slop.force.replace(x, slop.interval[1] - slop.interval[0]),
                    nextSlop.force.replace(x, 0))
            print('s{c}({d}) = s{cp}(0) --> {eq1} = {eq2}'.format(c=d_counter, eq1=eq.args[0], eq2=eq.args[1],
                                                                  cp=d_counter + 1,
                                                                  d=deflection.interval[1] - deflection.interval[
                                                                      0]))
            d_counter += 1
            equation.append(eq)

        solve = smp.solve(equation, cost)
        self.deflectionSlope = [
            Fun((1 / (self.I * self.E)) * i.force.subs(solve).replace(x, x + i.interval[0]), i.interval) for i in
            self.deflectionSlope]
        self.deflections = [Fun((1 / (self.I * self.E)) * i.force.subs(solve).replace(x, x + i.interval[0]), i.interval)
                            for i in self.deflections]

        # self.deflections = [ for i in self.deflections]
        for i in self.deflections:
            print(i.force)


class Fun:
    def __init__(self, load, interval):
        self.force = (smp.Abs(0) + load)
        self.interval = interval


class Moment:
    def __init__(self, force, offset):
        self.force = force
        self.offset = offset


class Reaction:
    def __init__(self, force, offset):
        self.force = force
        self.offset = offset


class Support:
    def __init__(self, _type, offset):
        self.type = _type
        self.offset = offset
