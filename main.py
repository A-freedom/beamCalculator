from beam import *

if __name__ == '__main__':
    _name = "Singer & Andre Problem 419"
    _length = 9
    _load = [Load((270 / 6) * x, [0, 6])]
    _reaction = [Reaction(450, 0), Reaction(360, 9)]
    _moments = []
    beam = Beam(_load, _reaction, _moments, _length, name=_name)
    # beam.printDetails()
    beam.displayPlots()

    user_input = input('do you want to save the beams to the test subjects?\ntype "save it" ::')
    if user_input == "save it":
        beam.save()
        print("the beam is saved 🥳")
    print('program is finished ,thanks 😘')
