from beam import *


def run_all():
    with open("lib/backup.txt", "r") as read_file:
        testing_data = json.load(read_file)
    for expect in testing_data['beams']:
        name = expect['name']
        reactions = [Reaction(re[0], re[1]) for re in expect['reactions']]
        moments = [Moment(me[0], me[1]) for me in expect['moments']]
        loads = [Load(smp.parse_expr(lo[0]), lo[1]) for lo in expect['loads']]
        l = expect['l']
        beam = Beam(loads, reactions, moments, l, name=name)
        beam.displayPlots()


if __name__ == '__main__':
    run_all()
    # _name = "beam number 1"
    # _length = 4
    # _load = [Load(7.5 * x, [0, 2]), Load(15 - 7.5 * x, [2, 4])]
    # _reaction = [Reaction(15, 1), Reaction(15, 3)]
    # _moments = []
    # beam = Beam(_load, _reaction, _moments, _length, name=_name)
    # beam.displayPlots()
    # beam.printDetails()
    #
    # user_input = input('do you want to save the beams to the test subjects?\ntype "save it" ::')
    # if user_input == "save it":
    #     beam.save()
    #     print("the beam is saved 🥳")
    # print('program is finished ,thanks 😘')
