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
        reference = 0.0
        heterozygous = 0.5
        homozygous = 1.0
        genotypes = [reference, heterozygous, homozygous]
        prior = 1.0 / 3.0
        nref = self.nref
        nalt = self.nalt
        pdg = {}
        for g in genotypes:
            # data likelihood P(D\g)
            pbin = scipy.stats.binom_test(nalt, nref + nalt, g, alternative='two-sided')
            pdg[g] = pbin
        posterior = {g: prior * pbin for g, pbin in pdg.items()}  # That's not actually the posterior, it's an intermediate before regularization
        regularization = sum([pbinp for pbinp in posterior.values()])
        posterior = {g: p / regularization for g, p in posterior.items()}
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
