
from decimal import Decimal

class peak(object):
    def __init__(self, mass, inten, charge):
        self.mass = mass
        self.inten = inten
        self.charge = charge
    def to_string(self):
        print self.mass, self.inten, self.charge

class spec(object):
    def __init__(self):
        self.id = 0
        self.scan = 0
        self.activation = ""
        self.pre_mz = Decimal(0.0)
        self.pre_charge = 0
        self.pre_mass = Decimal(0.0)
        self.peak_lst = []

    def add_scan(self, scan):
        self.scan = scan

    def add_act(self, act):
        self.activation = act

    def add_id(self, id):
        self.id = id

    def add_mz(self, mz):
        self.pre_mz = mz

    def add_charge(self, charge):
        self.pre_charge = charge

    def add_mass(self, mass):
        self.pre_mass = mass

    def add_peak(self, p):
        self.peak_lst.append(p)

    def to_file(self, filename):
        out = open(filename, "a")
        out.write("BEGIN IONS\n")
        out.write("ID=" + str(self.id) + "\n")
        out.write("SCANS=" + str(self.scan) + "\n")
        out.write("ACTIVATION=" + self.activation + "\n")
        out.write("PRECURSOR_MZ=" + str(self.pre_mz) + "\n")
        out.write("PRECURSOR_CHARGE=" + str(self.pre_charge) + "\n")
        out.write("PRECURSOR_MASS=" + str(self.pre_mass) + "\n")
        for p in self.peak_lst:
            out.write(str(p.mass) + "\t" + str(p.inten) + "\t" + str(p.charge) + "\n")
        out.write("END IONS\n\n")


def readMsAlign(filename):
    msalign_file = open(filename, 'r').readlines()

    idx = 0

    spec_lst = []

    while (idx < len(msalign_file)):
        line = msalign_file[idx].rstrip('\r\n')
        if line == "BEGIN IONS" :
            sp = spec()
            idx = idx + 1
            line = msalign_file[idx].rstrip('\r\n')
            sp.add_id(int(float(line.split("=")[1])))
            idx = idx + 1
            line = msalign_file[idx].rstrip('\r\n')
            sp.add_scan(int(float(line.split("=")[1])))
            idx = idx + 1
            line = msalign_file[idx].rstrip('\r\n')
            sp.add_act(line.split("=")[1])
            idx = idx + 1
            line = msalign_file[idx].rstrip('\r\n')
            sp.add_mz(Decimal(line.split("=")[1]))
            idx = idx + 1
            line = msalign_file[idx].rstrip('\r\n')
            sp.add_charge(int(float(line.split("=")[1])))
            idx = idx + 1
            line = msalign_file[idx].rstrip('\r\n')
            sp.add_mass(Decimal(line.split("=")[1]))
            idx = idx + 1
            while (idx < len(msalign_file)):
                line = msalign_file[idx].rstrip('\r\n')
                if line == "END IONS":
                    spec_lst.append(sp)
                    break;
                tmp = line.split("\t")
                p = peak(Decimal(tmp[0]), Decimal(tmp[1]), tmp[2])
                sp.add_peak(p)
                idx = idx + 1

        idx = idx + 1
    
    return spec_lst;
