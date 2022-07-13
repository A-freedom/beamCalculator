from beam import *

if __name__ == '__main__':
    # with open("beam_test.json", "r") as read_file:
    #     testing_data = json.load(read_file)
    # for expect in testing_data['beams']:
    _name = 'expect[]'
    _supports = [Support('pin', 0), Support('pin', 6)]
    _reactions = [Reaction(50, 0), Reaction(50, 6),Reaction(-20,2)]
    _moments = []
    _loads = [Fun(10, [0, 2]), Fun(15, [2, 6])]
    _length = 6
    _beam = Beam(_loads, _reactions, _moments, _length, name=_name, supports=_supports, E=1, I=1)
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    _beam.printDetails()
    # _beam.displayPlots()

    # user_input = input('do you want to save the beams to the test subjects?\ntype "save it" ::')
    # if user_input == "save it":
    #     beam.save()
    #     print("the beam is saved 🥳")

    print('program is finished ,thanks 😘')
