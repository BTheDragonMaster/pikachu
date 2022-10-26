from pikachu.errors import ColourError

BLACK = '#000000'
WHITE = '#ffffff'

LIGHT_GREY = "#CFCFCF"
GREY = "#B2B2B2"
DARK_GREY = "#727272"

LIGHT_BLUE = "#8FC3FA"
BLUE = "#3989DE"
DARK_BLUE = '#23268A'

LIGHT_RED = "#FF7E7E"
RED = "#DC3D3D"
DARK_RED = '#9B2727'

LIGHT_GREEN = "#63DD63"
GREEN = '#198719'
DARK_GREEN = "#125F12"

PINK = '#E766D0'
HOT_PINK = "#D722A6"
ORANGE = "#FF7400"
YELLOW = "#FFD100"
CYAN = "#00CEDB"
TURQUOISE = "#00AB9E"
PURPLE = "#61089E"
LIME = "#00E126"
RASPBERRY = "#E1005F"
BROWN = "#814724"

RANDOM_PALETTE_1 = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
                    '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
                    '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
                    '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080']

RANDOM_PALETTE_2 = [BLUE, RED, GREEN, PINK, ORANGE, CYAN, PURPLE, LIME, RASPBERRY,
                    DARK_RED, LIGHT_GREEN, DARK_BLUE, YELLOW, HOT_PINK, LIGHT_BLUE,
                    TURQUOISE, LIGHT_RED, DARK_GREEN]

string_to_colour = {"Red": RED,
                    "red": RED,
                    "Blue": BLUE,
                    "blue": BLUE,
                    "Green": GREEN,
                    "green": GREEN,
                    "Grey": GREY,
                    "Gray": GREY,
                    "grey": GREY,
                    "gray": GREY,
                    "Black": BLACK,
                    "black": BLACK,
                    "White": WHITE,
                    "white": WHITE,
                    "Orange": ORANGE,
                    "orange": ORANGE,
                    "Brown": BROWN,
                    "brown": BROWN,
                    "Cyan": CYAN,
                    "cyan": CYAN,
                    "Raspberry": RASPBERRY,
                    "raspberry": RASPBERRY,
                    "Purple": PURPLE,
                    "purple": PURPLE,
                    "Pink": PINK,
                    "pink": PINK,
                    "Lime": LIME,
                    "lime": LIME,
                    "Turquoise": TURQUOISE,
                    "turquoise": TURQUOISE,
                    "Yellow": YELLOW,
                    "yellow": YELLOW,
                    "hotpink": HOT_PINK,
                    "hot-pink": HOT_PINK,
                    "hot pink": HOT_PINK,
                    "Hot Pink": HOT_PINK,
                    "hotPink": HOT_PINK,
                    "Hot pink": HOT_PINK,
                    "hot_pink": HOT_PINK,
                    "lightred": LIGHT_RED,
                    "light-red": LIGHT_RED,
                    "light red": LIGHT_RED,
                    "Light Red": LIGHT_RED,
                    "lightRed": LIGHT_RED,
                    "Light red": LIGHT_RED,
                    "light_red": LIGHT_RED,
                    "lightblue": LIGHT_BLUE,
                    "light-blue": LIGHT_BLUE,
                    "light blue": LIGHT_BLUE,
                    "Light Blue": LIGHT_BLUE,
                    "lightBlue": LIGHT_BLUE,
                    "Light blue": LIGHT_BLUE,
                    "light_blue": LIGHT_BLUE,
                    "lightgreen": LIGHT_GREEN,
                    "light-green": LIGHT_GREEN,
                    "light green": LIGHT_GREEN,
                    "Light Green": LIGHT_GREEN,
                    "lightGreen": LIGHT_GREEN,
                    "Light green": LIGHT_GREEN,
                    "light_green": LIGHT_GREEN,
                    "darkred": DARK_RED,
                    "dark-red": DARK_RED,
                    "dark red": DARK_RED,
                    "Dark Red": DARK_RED,
                    "darkRed": DARK_RED,
                    "Dark red": DARK_RED,
                    "dark_red": DARK_RED,
                    "darkblue": DARK_BLUE,
                    "dark-blue": DARK_BLUE,
                    "dark blue": DARK_BLUE,
                    "Dark Blue": DARK_BLUE,
                    "darkBlue": DARK_BLUE,
                    "Dark blue": DARK_BLUE,
                    "dark_blue": DARK_BLUE,
                    "darkgreen": DARK_GREEN,
                    "dark-green": DARK_GREEN,
                    "dark green": DARK_GREEN,
                    "Dark Green": DARK_GREEN,
                    "darkGreen": DARK_GREEN,
                    "Dark green": DARK_GREEN,
                    "dark_green": DARK_GREEN,
                    "lightgrey": LIGHT_GREY,
                    "light-grey": LIGHT_GREY,
                    "light grey": LIGHT_GREY,
                    "Light Grey": LIGHT_GREY,
                    "lightGrey": LIGHT_GREY,
                    "Light grey": LIGHT_GREY,
                    "light_grey": LIGHT_GREY,
                    "lightgray": LIGHT_GREY,
                    "light-gray": LIGHT_GREY,
                    "light gray": LIGHT_GREY,
                    "Light Gray": LIGHT_GREY,
                    "lightGray": LIGHT_GREY,
                    "Light gray": LIGHT_GREY,
                    "light_gray": LIGHT_GREY,
                    "darkgrey": DARK_GREY,
                    "dark-grey": DARK_GREY,
                    "dark grey": DARK_GREY,
                    "Dark Grey": DARK_GREY,
                    "darkGrey": DARK_GREY,
                    "Dark grey": DARK_GREY,
                    "dark_grey": DARK_GREY,
                    "darkgray": DARK_GREY,
                    "dark-gray": DARK_GREY,
                    "dark gray": DARK_GREY,
                    "Dark Gray": DARK_GREY,
                    "darkGray": DARK_GREY,
                    "Dark gray": DARK_GREY,
                    "dark_gray": DARK_GREY,
                    }

def get_hex(colour):
    if colour.startswith('#') and len(colour) == 7:
        return colour
    elif colour in string_to_colour:
        return string_to_colour[colour]
    else:
        raise ColourError(colour)
