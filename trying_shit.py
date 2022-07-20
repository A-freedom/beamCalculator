from sympy import diff, solve
from sympy.abc import x


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


def findMixAndMin(function, interval):
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


d = findMixAndMin(
    -4.69219219219219e-10 * x ** 4 + 1.25125125125125e-9 * x ** 3 + 3.003003003003E-8 * x ** 2 - 5.92258925592259e-8 * x - 2.03536870203537e-7,
    [2, 6])
print(d)
