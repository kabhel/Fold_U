"""
.. module:: classes
   :synopsis: This module implements Alignment, Template and Residue classes
"""

# Third-party modules
import numpy as np



class Alignment:
    """
    .. class:: Alignment

      This class groups informations about an alignment.

    Attributes:
        score: Score of the alignment
        query_residues: Query's sequence of residues as list of Residues objects
        template: Instance of a Template object as
                  ``Template(template_name, template_residues)``
    """

    def __init__(self, score, query_residues, template_name, template_residues):
        """The constructor of an instance of the Alignment class."""
        self.score = score
        self.query_residues = query_residues
        self.template = Template(template_name, template_residues)


class Template:
    """
    .. class:: Template

      This class groups informations about a template sequence/structure.

    Attributes:
        name: Name of the template
        residues: Template's sequence of residues as list of Residues objects
        pdb: PDB code of the template
    """

    def __init__(self, name, residues):
        """The constructor of an instance of the Template class."""
        self.name = name
        self.residues = residues
        self.pdb = None

    def get_pdb(self, metafold_dict):
        """
            Get the PDB file name of the current template from the template's name.

            Args:
                self: The current alignment's template.
                dictionary: A dictionary with key = template name and value = pdb file

            Returns:
                void
        """
        self.pdb = metafold_dict[self.name]

    def get_all_ca_coords(self):
        """
            Parse the PDB file and sets the coordinates of the residues list of
            the current template.

            Args:
                self: The current alignment's template.

            Returns:
                void
        """
        res_num = 0
        with open("data/pdb/" + self.name + "/" + self.pdb, 'r') as file:
            for line in file:
                line_type = line[0:6].strip()
                name_at = line[12:16].strip()
                if line_type == "ATOM" and name_at == "CA":
                    x_coord = float(line[30:38].strip())
                    y_coord = float(line[38:46].strip())
                    z_coord = float(line[46:54].strip())
                    # Skip gaps in the template
                    while self.residues[res_num].name == '-':
                        res_num += 1
                    self.residues[res_num].ca_coords = np.array([x_coord, y_coord, z_coord])
                    res_num += 1


class Residue:
    """
    .. class:: Residue

      This class groups informations about a residue.

    Attributes:
        name: Name of the residue (1 letter code)
        ca_coords: 3D coordinates of the residue
    """

    def __init__(self, name):
        """The constructor of an instance of the Residue class."""
        self.name = name
        self.ca_coords = None

    def __str__(self):
        if self.ca_coords is None:
            return "<" + self.name + " | " + "empty coordinates>"
        coords = [str(coord) for coord in self.ca_coords]
        return "<" + self.name + " | " + str(coords) + ">"