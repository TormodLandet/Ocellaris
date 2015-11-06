class OcellarisError(Exception):
    def __init__(self, header, description):
        super(OcellarisError, self).__init__('%s: %s' % (header, description))
        self.header = header
        self.description = description

def report_error(header, description):
    raise OcellarisError(header, description)
    