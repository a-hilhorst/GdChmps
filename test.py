import numpy as np
import matplotlib.pyplot as plt
import PySimpleGUI as Sg

layout = [
    [Sg.Text('# images of each series: ',
             tooltip='Enter the number of images of each series as an integer separated by a comma \n'
                     'e.g. 2,3,4 for 3 series of 2, 3, and 4 images respectively.',
             size=(20, 1)), Sg.Input(key='n_im')],
    [Sg.Submit('T-E-S-T'), Sg.Exit()]
]

window = Sg.Window('Test', layout)

while True:
    event, values = window.read()
    if event in (None, 'Exit'):
        break

    n_im = []
    while True:
        try:
            n_im = values['n_im'].split(',')
            n_im = [i.split('x') for i in n_im]
            n_im = [[int(i) for i in j] for j in n_im]
        except ValueError:
            Sg.Popup('Error', 'Please enter # of images as integers separated by a comma\n'
                              'Please enter magnification as doubles separated by a comma')
            event, values = window.read()
            continue
        break

    if event in 'T-E-S-T':
        Sg.Popup(n_im)

