from beam import *


def run_all():
    with open("lib/backup.txt", "r") as read_file:
        testing_data = json.load(read_file)
    for expect in testing_data['beams']:
        name = expect['name']
        reactions = [Reaction(re[0], re[1]) for re in expect['reactions']]
        moments = [Moment(me[0], me[1]) for me in expect['moments']]
        loads = [Fun(smp.parse_expr(lo[0]), lo[1]) for lo in expect['loads']]
        l = expect['l']
        beam = Beam(loads, reactions, moments, l, name=name)
        beam.displayPlots()


if __name__ == '__main__':
    _name = "deflection ex 6"
    _length = 5
    _supports = [Support('pin', 1), Support('pin', 4)]
    _load = [Fun(12, [1, 4])]
    _reaction = [Reaction(18, 1), Reaction(18, 4)]
    _moments = []
    beam = Beam(_load, _reaction, _moments, _length, name=_name, supports=_supports, E=200E9, I=3500E-8)
    beam.printDetails()

    # print('the maximum deflection::',beam.deflection[0].force.evalf(subs={x:1.5}))
    beam.displayPlots()

    # user_input = input('do you want to save the beams to the test subjects?\ntype "save it" ::')
    # if user_input == "save it":
    #     beam.save()
    #     print("the beam is saved 🥳")

    print('program is finished ,thanks 😘')
