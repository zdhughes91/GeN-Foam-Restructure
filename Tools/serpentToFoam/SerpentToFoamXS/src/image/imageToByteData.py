"""
SerpentToFoamXS

Author: Thomas Guilbaud, EPFL/Transmutex SA
Last Update: 18/03/2022
"""

import base64
import os

# Translate image into byte data
def pic2str(file, functionName):
    picture = open(file, 'rb')
    content = '{} = {}\n'.format(functionName, base64.b64encode(picture.read()))
    picture.close()

    with open('./SerpentToFoamXS/image/logoB64.py', 'w') as f:
        f.write(content)

if __name__ == '__main__':
    pic2str("SerpentToFoamXS/image/logo.png", "imageByte")
