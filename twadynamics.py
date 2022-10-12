import numpy as np
from numpy import linalg as LA


class evo:
    def __init__(self, cavityfreq=1600, vibrationfreq=1600, huangrhys=0.75, jccoupling=0.01e-3):
        wavenumberToHartree = 4.55633e-6
        self.freqC = cavityfreq * wavenumberToHartree
        self.freqV = vibrationfreq * wavenumberToHartree
        self.hr = huangrhys
        self.g = jccoupling
        self.fieldEnDiag = np.array([])
        self.fieldSize = 0
        self.bathEnDiag = np.array([])
        self.bathSize = 0
        self.systemSize = 2
        self.eigvals = []
        self.eigvecs = []
        self.dfrequency = 0
        self.bCouplingFn = []
        self.fCouplingFn = []
        self.totalSize = self.systemSize
        self.hessian = np.zeros((self.totalSize, self.totalSize))

    def addbath(self, bathsize=500, fieldsize=500, dfreq=0.004, bdamp=0, fdamp=0):
        # add the bath and fields to the total size
        self.dfrequency = dfreq
        self.totalSize = self.systemSize + bathsize + fieldsize
        self.hessian = np.zeros((self.totalSize, self.totalSize))
        self.fieldSize = fieldsize
        self.bathSize = bathsize
        # Diagonal elements of the free field energies
        if fieldsize != 0:
            self.fieldEnDiag = np.linspace(self.freqC - dfreq / 2, self.freqC + dfreq / 2, fieldsize) ** 2
        # and bath energies
        if bathsize != 0:
            self.bathEnDiag = np.linspace(self.freqV - dfreq / 2, self.freqV + dfreq / 2, bathsize) ** 2

        # create the frequency dependent coupling function for both the external free fields
        self.fCouplingFn = np.multiply(2 * np.sqrt(self.freqC * np.sqrt(self.fieldEnDiag)),
                                       np.sqrt(fdamp * np.sqrt(self.fieldEnDiag) * dfreq / fieldsize))
        # self.fCouplingFn = np.full(self.fieldSize, 2 * np.sqrt(self.freqC * self.freqC) * (self.freqC/(2 * np.pi * (500/0.004) * 1000000))**0.5)
        # and bath
        self.bCouplingFn = np.multiply(2 * np.sqrt(self.freqV * np.sqrt(self.bathEnDiag)),
                                       np.sqrt(bdamp * np.sqrt(self.bathEnDiag) * dfreq / bathsize))

    def updatefCoupling(self, fdamp):
        self.fCouplingFn = np.multiply(2 * np.sqrt(self.freqC * np.sqrt(self.fieldEnDiag)),
                                       np.sqrt(fdamp * np.sqrt(self.fieldEnDiag) * self.dfrequency / self.fieldSize))

    def updatebCoupling(self, bdamp):
        self.bCouplingFn = np.multiply(2 * np.sqrt(self.freqV * np.sqrt(self.bathEnDiag)),
                                       np.sqrt(bdamp * np.sqrt(self.bathEnDiag) * self.dfrequency / self.bathSize))

    def occ(self, time, cav=True, vib=True, field=True, bath=True):
        # Determine the re-organization energy, lambda.
        lmbda = self.hr * self.freqV
        # Calculates the Hessain matrix coupling strength from the Jaynes-Cummings coupling strength.
        G = 2 * np.sqrt(self.freqC * self.freqV) * self.g
        # Build 2x2 system self.hessian matrix
        self.hessian[0, 0] = self.freqC ** 2
        self.hessian[1, 1] = self.freqV ** 2
        self.hessian[1, 0] = G
        self.hessian[0, 1] = G

        # bath component of self.hessian matrix
        i = 0
        for en in self.fieldEnDiag:
            self.hessian[self.systemSize + i, self.systemSize + i] = en
            self.hessian[0, self.systemSize + i] = self.fCouplingFn[i]
            self.hessian[self.systemSize + i, 0] = self.fCouplingFn[i]
            i += 1

        i = 0
        for en in self.bathEnDiag:
            self.hessian[self.systemSize + self.fieldSize + i, self.systemSize + self.fieldSize + i] = en
            self.hessian[1, self.systemSize + self.fieldSize + i] = self.bCouplingFn[i]
            self.hessian[self.systemSize + self.fieldSize + i, 1] = self.bCouplingFn[i]
            i += 1

        # Calculate the eigenvalues and eigenvectors for the self.hessian. u is the unitary matrix that diagonalizes the
        # self.hessian.
        hvals, u = LA.eigh(self.hessian)
        self.eigvals = hvals
        self.eigvecs = u

        # calculated constant prefactors

        # cavity:
        if cav:
            acav1 = 1 / (2 * self.freqC)
            acav2 = 2 * self.freqC
            acav3 = 1 / self.freqC
            acav4 = self.freqC
        # vibration:
        if vib:
            avib1 = 1 / (2 * self.freqV)
            avib2 = 2 * self.freqV
            avib3 = 1 / self.freqV
            avib4 = self.freqV
        # fields
        afield1 = 1 / (2 * np.sqrt(self.fieldEnDiag))
        afield2 = 2 * np.sqrt(self.fieldEnDiag)
        afield3 = 1 / np.sqrt(self.fieldEnDiag)
        afield4 = np.sqrt(self.fieldEnDiag)
        abath1 = 1 / (2 * np.sqrt(self.bathEnDiag))
        abath2 = 2 * np.sqrt(self.bathEnDiag)
        abath3 = 1 / np.sqrt(self.bathEnDiag)
        abath4 = np.sqrt(self.bathEnDiag)

        # cavity:
        if cav:
            # define time dependent matrices
            ccav1 = np.zeros((hvals.size, time.size))
            ccav2 = np.zeros((hvals.size, time.size))

            # define matrices for the thermal term.
            ccav3 = np.sqrt(hvals)
            ccav4 = 1 / ccav3

            # transformation coefficients for transforming from localized basis to delocalized polariton basis
            bcav = u[0, :]

        if vib:
            cvib1 = np.zeros((hvals.size, time.size))
            cvib2 = np.zeros((hvals.size, time.size))

            cvib3 = np.sqrt(hvals)
            cvib4 = 1 / cvib3

            bvib = u[1, :]

        if field:
            cfield1 = np.zeros((hvals.size, time.size))
            cfield2 = np.zeros((hvals.size, time.size))

            cfield3 = np.sqrt(hvals)
            cfield4 = 1 / cfield3

            bfield = u[self.systemSize:self.fieldSize + self.systemSize, :]

        if bath:
            cbath1 = np.zeros((hvals.size, time.size))
            cbath2 = np.zeros((hvals.size, time.size))

            cbath3 = np.sqrt(hvals)
            cbath4 = 1 / cbath3

            bbath = u[self.systemSize + self.fieldSize:self.systemSize + self.fieldSize + self.bathSize, :]


        # calculation of occupation values as a function of time using matrix math and numpy for efficiency
        i = 0
        for Omega in np.sqrt(hvals):
            # multiply the current eigen frequency by each point in the time array
            omeg_t = Omega * time

            if cav:
                ccav1[i] = u[1, i] * np.sin(omeg_t) / Omega
                ccav2[i] = u[1, i] * (np.sin(omeg_t / 2) ** 2) / Omega ** 2

            if vib:
                cvib1[i] = u[1, i] * np.sin(omeg_t) / Omega
                cvib2[i] = u[1, i] * (np.sin(omeg_t / 2) ** 2) / Omega ** 2

            if field:
                cfield1[i] = u[1, i] * np.sin(omeg_t) / Omega
                cfield2[i] = u[1, i] * (np.sin(omeg_t / 2) ** 2) / Omega ** 2

            if bath:
                cbath1[i] = u[1, i] * np.sin(omeg_t) / Omega
                cbath2[i] = u[1, i] * (np.sin(omeg_t / 2) ** 2) / Omega ** 2
            i += 1

        # matrix math that calculates occupation values as a function of time
        ncav = 0
        nvib = 0
        nfield = 0
        nbath = 0
        if cav:
            ncav1 = acav1 * np.square(bcav.dot(ccav1)) + acav2 * np.square(bcav.dot(ccav2))
            ncav2 = acav3 * np.square(bcav).dot(ccav3) + acav4 * np.square(bcav).dot(ccav4)
            ncav = 2 * lmbda * self.freqV ** 2 * ncav1 + (1 / 4) * ncav2 - 0.5
        if vib:
            nvib1 = avib1 * np.square(bvib.dot(cvib1)) + avib2 * np.square(bvib.dot(cvib2))
            nvib2 = - 4 * lmbda * self.freqV * bvib.dot(cvib2)
            nvib3 = avib3 * np.square(bvib).dot(cvib3) + avib4 * np.square(bvib).dot(cvib4)
            nvib = 2 * lmbda * self.freqV ** 2 * nvib1 + nvib2 + (1 / 4) * nvib3 - 0.5 + (lmbda / self.freqV)
        if field:
            nfield1 = afield1.dot(np.square(bfield.dot(cfield1))) + afield2.dot(np.square(bfield.dot(cfield2)))
            nfield2 = afield3.dot(np.square(bfield).dot(cfield3)) + afield4.dot(np.square(bfield).dot(cfield4))
            nfield = 2 * lmbda * self.freqV ** 2 * nfield1 + (1 / 4) * nfield2 - self.fieldSize * 0.5
        if bath:
            nbath1 = abath1.dot(np.square(bbath.dot(cbath1))) + abath2.dot(np.square(bbath.dot(cbath2)))
            nbath2 = abath3.dot(np.square(bbath).dot(cbath3)) + abath4.dot(np.square(bbath).dot(cbath4))
            nbath = 2 * lmbda * self.freqV ** 2 * nbath1 + (1 / 4) * nbath2 - self.bathSize * 0.5
        return ncav, nvib, nfield, nbath

    def spectra(self, time):
        # Determine the re-organization energy, lambda.
        lmbda = self.hr * self.freqV
        # Calculates the Hessain matrix coupling strength from the Jaynes-Cummings coupling strength.
        G = 2 * np.sqrt(self.freqC * self.freqV) * self.g
        # Build 2x2 system self.hessian matrix
        self.hessian[0, 0] = self.freqC ** 2
        self.hessian[1, 1] = self.freqV ** 2
        self.hessian[1, 0] = G
        self.hessian[0, 1] = G

        # bath component of self.hessian matrix
        envCoupling = np.zeros((self.totalSize, self.totalSize))
        i = 0
        for en in self.fieldEnDiag:
            self.hessian[self.systemSize + i, self.systemSize + i] = en
            envCoupling[0, self.systemSize + i] = self.fCouplingFn[i]
            i += 1

        i = 0
        for en in self.bathEnDiag:
            self.hessian[self.systemSize + self.fieldSize + i, self.systemSize + self.fieldSize + i] = en
            envCoupling[1, self.systemSize + self.fieldSize + i] = self.bCouplingFn[i]
            i += 1

        self.hessian = self.hessian + envCoupling + envCoupling.T

        # Calculate the eigenvalues and eigenvectors for the self.hessian. u is the unitary matrix that diagonalizes the
        # self.hessian.
        hvals, u = LA.eigh(self.hessian)
        self.eigvals = hvals
        self.eigvecs = u

        # calculated constant prefactors
        afield1 = np.nan_to_num(1 / (2 * np.sqrt(self.fieldEnDiag)), posinf=0)
        afield2 = np.nan_to_num(2 * np.sqrt(self.fieldEnDiag), posinf=0)
        afield3 = np.nan_to_num(1 / np.sqrt(self.fieldEnDiag), posinf=0)
        afield4 = np.sqrt(self.fieldEnDiag)

        # define time dependent matrices
        cfield1 = np.zeros((hvals.size, time.size))
        cfield2 = np.zeros((hvals.size, time.size))

        # define matrices for the thermal term.
        cfield3 = np.sqrt(hvals)
        cfield4 = 1 / cfield3

        # calculation of occupation values as a function of time using matrix math and numpy for efficiency
        i = 0
        for Omega in np.sqrt(hvals):
            omeg_t = Omega * time
            cfield1[i] = u[1, i] * np.sin(omeg_t) / Omega
            cfield2[i] = u[1, i] * (np.sin(omeg_t / 2) ** 2) / Omega ** 2
            i += 1

        spectra = np.zeros(self.fieldEnDiag.size)
        for i in range(self.fieldEnDiag.size):
            # transformation coefficients for transforming from localized basis to delocalized polariton basis
            bfield = u[self.systemSize + i, :]

            # matrix math that calculates occupation values as a function of time
            nfield1 = afield1[i] * np.square(bfield.dot(cfield1)) + afield2[i] * np.square(bfield.dot(cfield2))
            nfield2 = afield3[i] * np.square(bfield).dot(cfield3) + afield4[i] * np.square(bfield).dot(cfield4)
            nfield = 2 * lmbda * self.freqV ** 2 * nfield1 + (1 / 4) * nfield2 - 0.5
            spectra[i] = nfield[-1]

        return spectra
