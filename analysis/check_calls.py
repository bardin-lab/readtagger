"""Goes through a list of images. The user can press 'n' or 'y',
which will cause True or False to be written in the output file as a column next to the filename.
Prerequisite is a function opencv library (on OSX can be installed by brew install opencv), pyobjc (pip install pyobjc)
and opencv-python
"""

import os
import click
import cv2
from AppKit import NSScreen


def check_image(filename, height):
    image = cv2.imread(filename)
    image_name = os.path.basename(filename)
    is_positive = _check(image, image_name, height)
    return image_name, is_positive


def _check(image, image_name, target_height):
    cv2.namedWindow(image_name, cv2.WINDOW_NORMAL)
    while True:
        height, width = image.shape[:2]
        scaling_f = height / target_height
        small = cv2.resize(image, None, fx=1 / scaling_f, fy=1 / scaling_f)
        # display the image and wait for a keypress
        cv2.imshow(image_name, small)
        key = cv2.waitKey(1)

        # if the 'y' key is pressed return True
        if key == ord("y"):
            cv2.destroyWindow(image_name)
            print('True')
            return True

        # if the 'n' key is pressed return False
        elif key == ord("n"):
            print('False')
            cv2.destroyWindow(image_name)
            return False

        # if the 'm' key is pressed return Maybe
        elif key == ord("m"):
            print('Maybe')
            cv2.destroyWindow(image_name)
            return 'Maybe'


@click.command()
@click.argument('files',
                nargs=-1, type=click.Path(exists=True))
@click.option('--output_path',
              help='Write result here',
              default=None,
              type=click.Path(exists=False))
def check_files(files, output_path):
    height = NSScreen.mainScreen().frame().size.height
    with open(output_path, 'w') as out:
        for file in files:
            image_name, is_positive = check_image(file, height)
            template = "%s\t%s\n"
            out.write(template % (image_name, is_positive))


if __name__ == "__main__":
    check_files()
