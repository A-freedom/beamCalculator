from sympy import SingularityFunction, Abs, integrate
from sympy.abc import x


class Singularity:
    def __init__(self, ex, a):
        self.a = a
        self.ex = ex

    def evalf(self, X):
        if X - self.a >= 0:
            return self.ex.evalf(subs={x: X - self.a})
        return 0

    def getIntegrate(self):
        return Singularity(integrate(self.ex, x), self.a)


if __name__ == '__main__':
    f = Singularity(10 - x, 2)
    for i in range(-5, 5):
        print(i, '-->', f.evalf(i))
