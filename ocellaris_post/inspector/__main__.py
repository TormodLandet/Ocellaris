import sys, os
from .wx_gui import show_inspector


def run_from_console():
    """
    Parse command line arguments and run the wxPython GUI
    """
    # Get report files to save
    h5_file_names = sys.argv[1:]
    
    # Get lables
    lables = []
    for i in range(len(h5_file_names)):
        fn = h5_file_names[i]
        if ':' in fn:
            fn, label = fn.split(':')
            h5_file_names[i] = fn
        else:
            bname = os.path.basename(fn)
            bname_split = bname.split('_endpoint_')
            label = bname_split[0]
        lables.append(label)
    
    # Make plots
    show_inspector(h5_file_names, lables)


if __name__ == '__main__':
    run_from_console()
