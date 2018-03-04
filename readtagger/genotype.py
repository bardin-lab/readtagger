"""Genotype module."""
import sys
import scipy.stats


class Genotype(list):
    """A Genotype object."""

    def __init__(self, nref, nalt):
        """
        Genotype object.

        nref is the amount of evidence supporting the reference allele,
        nalt is the evidence supporting an alternative allele.

        >>> Genotype(13750, 5257).genotype
        'homozygous'
        """
        self.nref = nref
        self.nalt = nalt
        self.genotype_likelihood()

    @property
    def reference(self):
        """Return p-value decribing the probability that the genotype is reference."""
        if len(self) == 3:
            return self[0]

    @property
    def heterozygous(self):
        """Return p-value decribing the probability that the genotype is heterozygous."""
        if len(self) == 3:
            return self[1]

    @property
    def homozygous(self):
        """Return p-value decribing the probability that the genotype is homozygous."""
        if len(self) == 3:
            return self[2]

    def genotype_likelihood(self):
        r"""
        Calculate genotype likelihood.

        P(g|D) = P(g)P(D\g)/sum(P(g)P(D|g')) where P(D|g) = Pbin(Nalt, Nalt + Nfef)
        :return:
        """
        reference = 0.03
        heterozygous = 0.5
        homozygous = 0.97
        genotypes = [reference, heterozygous, homozygous]
        priors = [0.9, 0.05, 0.05]
        nref = self.nref
        nalt = self.nalt
        pdg = {}
        for g, prior in zip(genotypes, priors):
            # data likelihood P(D\g)
            pbin = scipy.stats.binom_test(nalt, nref + nalt, g, alternative='two-sided')
            pdg[g] = pbin * prior
        regularization = sum([pbinp for pbinp in pdg.values()])
        if regularization == 0:
            # This can happen if regularization is rounded to 0
            regularization += sys.float_info.min
        posterior = {g: p / regularization for g, p in pdg.items()}
        self.append(posterior[reference])
        self.append(posterior[heterozygous])
        self.append(posterior[homozygous])
        genotype_p = max([self.reference, self.heterozygous, self.homozygous])
        if genotype_p == self.homozygous:
            genotype = 'homozygous'
        elif genotype_p == self.heterozygous:
            genotype = 'heterozygous'
        elif genotype_p == self.reference:
            genotype = 'reference'
        self.genotype = genotype
