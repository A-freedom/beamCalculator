from beam import *

if __name__ == '__main__':
    # with open("beam_test.json", "r") as read_file:
    #     testing_data = json.load(read_file)
    # for expect in testing_data['beams']:
    _name = 'expect[]'
    _supports = [Support('pin', 1), Support('pin', 5)]
    _reactions = [PointedLoad(-20,2)]
    _moments = []
    _loads = [DistributedLoad(0, 2, 10, 10),DistributedLoad(2, 6, 15, 15)]
    _length = 6
    _beam = Beam(_loads, _reactions, _moments, _length, name=_name, supports=_supports, I=6660000, E=200)
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    _beam.printDetails()
    _beam.displayPlots()

    # user_input = input('do you want to save the beams to the test subjects?\ntype "save it" ::')
    # if user_input == "save it":
    #     beam.save()
    #     print("the beam is saved 🥳")

    print('program is finished ,thanks 😘')
