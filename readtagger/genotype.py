"""Genotype module."""
import logging
import scipy.stats


class Genotype(list):
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

    @property
    def reference(self):
        """Return p-value decribing the probability that the genotype is reference."""
        if len(self) == 3:
            return self[0]
        else:
            return None

    @property
    def heterozygous(self):
        """Return p-value decribing the probability that the genotype is heterozygous."""
        if len(self) == 3:
            return self[1]
        else:
            return None

    @property
    def homozygous(self):
        """Return p-value decribing the probability that the genotype is homozygous."""
        if len(self) == 3:
            return self[2]
        else:
            return None

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
        else:
            logging.info("Could not determine genotype, nref: %s, nalt: %s", self.nref, self.nalt)
            genotype = 'NA'
        self.genotype = genotype
