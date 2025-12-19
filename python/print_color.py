from enum import Enum

class Color(Enum):
    DEFAULT = 0
    Black = 16
    Red = 202
    Green = 46
    Yellow = 226
    White = 255
    Blue = 21
    Cyan = 51
    Magenta = 201


def print_color (text : str, c : Color = Color.DEFAULT, keep_color : bool = False, end_ : str='\n') :
    # on check si la couleur demandé est DEFAULT
    if c.value == 0 :
        print("\x1B[0m", end='')
    else :
        #préparation du motif pour écrire en couleur
        motif : str = f"\x1B[38;5;{c.value}m"
        print(motif, end='')  # affichage de la couleur :


    print(text, end=end_) # affichage du texte en console

    # on remet la console en mode normal si c'est demandé
    if not keep_color :
        print("\033[0m", end='')