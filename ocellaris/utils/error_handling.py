"""
Error handling and assertions giving nice user input
and controlled shutdown of Ocellaris
"""


class OcellarisError(Exception):
    def __init__(self, header, description):
        super(OcellarisError, self).__init__('%s: %s' % (header, description))
        self.header = header
        self.description = description


def ocellaris_error(header, description):
    raise OcellarisError(header, description)


def verify_key(name, key, options, loc=None):
    """
    Verify that a key is among a set of options. If not
    give a sensible warning.
    
    * name should be non-capitalized, ie. 'flower'
    * key should be the user provided input, ie. 'dandelion'
    * options should be allowable inputs, ie. ['rose', 'daisy']
    * loc is optional to provide more context, ie. 'get_flower2()'
    """
    if key not in options:
        loc = ' in %s' % loc if loc is not None else ''
        if len(options) > 1:
            if hasattr(options, 'keys'):
                options = list(options.keys())
            available_options = '\n'.join(' - %r' % m for m in options)
            ocellaris_error('Unsupported %s' % name,
                            'The %s %r is not available%s, please use one of:\n%s' %
                            (name, key, loc, available_options))
        else:
            available_options = ', '.join('%r' % m for m in options)
            ocellaris_error('Unsupported %s' % name,
                            'The %s %r is not available%s, only %s is available' %
                            (name, key, loc, available_options))
