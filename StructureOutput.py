from ase.io import write

class StructureOutput():
    """
    Class Output the Constructed Slab Structure with the help of ASE

    Attributes:
            @StructureWriter(): Output the constructed structure with a specific format
    """

    def __init__(self,filename,fileformat,structrue):

        self.filename = filename
        self.fileformat = fileformat
        self.structure = structrue

    def StructureWriter(self):

        write(self.filename, self.structure, format=self.fileformat)

        return