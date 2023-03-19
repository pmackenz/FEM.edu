"""
Defines an Example class from which all provided examples will be derived.

"""


class Example():

    def __init__(self):
        pass

    def __str__(self):
        s  = "Example {}\n".format(self.__class__.__name__)
        #s += len(s)*"-" + '\n'
        doc_string = self.docString()
        if doc_string:
            s += doc_string + '\n'
        return s

    def run(self):
        self.problem()
        return

        try:
            self.problem()
        except:
            print("failure to run Example: {}".format(self.__class__.__name__))

    def docString(self):
        """
        Return problem specific documentation as a multi-line string.
        """
        return ""

    def problem(self, *args, **kwargs):
        """
        Build and run the example problem
        """

        msg = """
        Missing problem() method for Example {}
        """.format(self.__class__.__name__)
        print(msg)


