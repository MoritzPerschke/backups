import requests
import json
from time import sleep

bridge = "http://192.168.0.116/api/4Qye7DLD4EY0jkFR7LnIPZiBF7ibJvhSmBYMFYOC"

with open('colors.json') as colorflie:
    colors = json.load(colorflie)


def change_color_single(light : int, color : dict):
    r = requests.put(bridge + f'/lights/{light}/state', data = json.dumps(color))
    if(not r.status_code == '200'):
        print(r.status_code)
        print(r.text)

def set_color(color : dict):
    for i in range(1,4):
        change_color_single(i, color)

def make_red():
    set_color(colors["red"])

def make_bright():
    set_color(colors["bright"])

while (True):
    make_red()
    sleep(0.1)
    make_bright()
    sleep(.1)