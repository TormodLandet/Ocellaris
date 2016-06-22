class OcellarisError(Exception):
    def __init__(self, header, description):
        super(OcellarisError, self).__init__('%s: %s' % (header, description))
        self.header = header
        self.description = description

def ocellaris_error(header, description):
    raise OcellarisError(header, description)

def verify_key(name, key, options, loc=None):
    if not key in options:
            available_otions = '\n'.join(' - %s' % m for m in options)
            loc = ' in %s' % loc if loc is not None else ''
            ocellaris_error('Unsupported %s' % name,
                            'The %s %r not available%s, please use one of:\n%s' %
                            (name, key, loc, available_otions))
