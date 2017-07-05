import sys
from . import show_inspector


def main():
    """
    Parse command line arguments and run the wxPython GUI
    """
    # Get report files to read
    h5_file_names = sys.argv[1:]
    
    # Get lables from command line arguments
    lables = []
    for i in range(len(h5_file_names)):
        fn = h5_file_names[i]
        if ':' in fn:
            fn, label = fn.split(':')
            h5_file_names[i] = fn
        else:
            label = None
        lables.append(label)
    
    # Show the Ocellaris Inspector
    show_inspector(h5_file_names, lables)


if __name__ == '__main__':
    main()
