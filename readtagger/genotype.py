"""Genotype module."""
import scipy.stats


class Genotype(object):
    """A Genotype object."""

    def __init__(self, nref, nalt):
        """
        Genotype object.

        nref is the amount of evidence supporting the reference allele,
        nalt is the evidence supporting an alternative allele.
        """
        self.nref = nref
        self.nalt = nalt
        self.genotype_likelihood()

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
        posterior = {g: p / regularization for g, p in pdg.items()}
        self.reference = posterior[reference]
        self.heterozygous = posterior[heterozygous]
        self.homozygous = posterior[homozygous]
        genotype_p = max([self.reference, self.heterozygous, self.homozygous])
        if genotype_p == self.homozygous:
            genotype = 'homozygous'
        if genotype_p == self.heterozygous:
            genotype = 'heterozygous'
        if genotype_p == self.reference:
            genotype = 'reference'
        self.genotype = genotype
