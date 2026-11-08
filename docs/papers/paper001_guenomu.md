# A Bayesian Supertree Model for Genome-Wide Species Tree Reconstruction

**Leonardo de Oliveira Martins, Diego Mallo and David Posada**
*Systematic Biology* **65**(3), 397–416 (2016).

[Published article and figures](https://doi.org/10.1093/sysbio/syu082) · [Software chapter](paper002_guenomu.md)

> Reading edition of the supplied manuscript. Scientific text and references are retained. Figures are available in the published article linked above; references to them open that article. The publisher’s version is authoritative.

## Abstract
Current phylogenomic data sets highlight the need for species tree methods able to deal with several sources of gene tree/species tree incongruence. At the same time, we need to make most use of all available data. Most species tree methods deal with single processes of phylogenetic discordance, namely, gene duplication and loss, incomplete lineage sorting or horizontal gene transfer. In this manuscript we address for the first time the problem of species tree inference from multilocus, genome-wide data sets in the presence of gene duplication and loss and incomplete lineage sorting, therefore without the need to identify orthologs or to use a single individual per species. We do this by extending the idea of Maximum Likelihood supertrees to a hierarchical Bayesian model where several sources of gene tree/species tree disagreement can be accounted for in a modular manner. We implemented this model in a computer program called *guenomu* whose input consists of posterior distributions of unrooted gene tree topologies for multiple gene families, and whose output is the posterior distribution of rooted species tree topologies. We conducted extensive simulations under complex phylogenomic models in order to evaluate the performance of our approach in comparison with other species tree approaches able to deal with multilabeled trees. Our method ranked best under simulated data sets, in spite of ignoring branch lengths, and performed well on empirical data sets as well. Our results show in addition that under complex simulation scenarios , gene tree parsimony is also a competitive approach once we consider its speed in contrast to more sophisticated models.

(Keywords: supertree, reconciliation, hierarchical Bayesian model, tree distance, phylogenomics)

## Introduction

The evolutionary history of a gene family is not an exact representation of the evolution of the species embedding it, since we must take into account the effects of gene duplication, population factors like incomplete lineage sorting as well as the possibility of horizontal gene transfer (Maddison 1997). In recent years the rapid accumulation of phylogenomic data has highlighted the need for specific methods to infer species phylogenies. Methods of species tree inference based on collection of gene alignments/trees can be broadly classified into supermatrix, supertree and model-based approaches. The supermatrix strategy consists in concatenating different genes into a single large alignment, assuming that minor conflicting phylogenetic signals will cancel out and a common phylogenetic trend -- the species history -- will emerge. Departures from this assumption must be handled by subjectively removing columns or whole regions (M. Wu, Chatterji, and Eisen 2012; Lanier, Huang, and Knowles 2014). This method relies on the correct inference of the orthologous sets, such that the sequences from a single species can be concatenated into a single 'supergene' (Degiorgio and Degnan 2010).

Alternatively, phylogenetic inference can be conducted independently for each gene alignment, and then in a second step the species tree is inferred from the resulting gene trees (Cotton and Wilkinson 2009). In this case the alignments are used indirectly, through summary statistics like maximum likelihood trees, to build a 'supertree'. Under the supertree umbrella, some methods try to find the species tree that minimizes its overall disagreement with the gene tree collection without explicitly taking into account the biological phenomena behind the incongruence, like the RF supertree (Mukul S. Bansal et al. 2010; Chaudhary, Burleigh, and Fernández-Baca 2013) and the MRP approaches (Ragan 1992; Bininda-Emonds 2004). Other supertree approaches aim for the most parsimonious scenario by assuming that the disagreement between the species tree and each of the gene trees is due to HGT (Whidden, Zeh, and Beiko 2013), duplications and losses or deep coalescences (Mukul S. Bansal and Eulenstein 2013). The problem of finding the species tree that minimizes the reconciliation cost is called Gene Tree Parsimony (GTP) (Mukul S. Bansal and Eulenstein 2013).

Most supertree approaches neglect the gene tree branch lengths, but there is another class of species tree methods closely related to the supertree methods which are based on the multispecies coalescent. These methods try to reconstruct the species tree from a matrix of distances between species, which in turn is built based on distance matrices from the individual gene trees (reviewed in (Liang Liu et al. 2009; Helmkamp, Jewett, and Rosenberg 2012)). These methods are based on the observation that the coalescence time of two genes from distinct species always precedes the speciation event. In general they assume that the gene trees are known without error (Kubatko, Carstens, and Knowles 2009; Liang Liu and Yu 2010). Some of these methods also can take into account hybridization (Meng and Kubatko 2009) or HGT (Yu et al. 2011).

In principle, the most powerful species tree methods are based on an explicit probabilistic modelling of the incongruence between species and gene trees based on gene alignments. Together with the phylogenetic likelihood that takes into account the stochastic uncertainty in the gene tree inference, these models also describe the probability of such gene tree being generated by a given species tree. If we assume that the sequences represent individuals from a given orthologous region, then this probability is given by the multispecies coalescent, that describes how lineages coalesce conditioned on a species tree (Rannala and Yang 2003). Under this model Bayesian methods have already been implemented to estimate the posterior distribution of species trees, together with ancestral population sizes and divergence times (Edwards, Liu, and Pearl 2007; Heled and Drummond 2010; Ronquist et al. 2012). If, however, we assume that each alignment represents a set of paralogous gene copies from a gene family, then we can model the duplications and losses inside each gene family through a birth-death process (Arvestad et al. 2003). Although theoretically it is possible to estimate the distribution of species trees while integrating over the particular birth-death histories, current implementations assume that the species tree is known (Akerborg et al. 2009; Rasmussen and Kellis 2011; Sjostrand et al. 2012) or try to find it by maximum likelihood (Boussau et al. 2013). Under the fixed tree assumption, it is furthermore possible to add HGT to the birth-death model (Szollosi et al. 2013). Fully probabilistic models, although more realistic, tend to be slow. Furthermore, they are usually limited to particular sources of gene tree / species tree disagreement. Under the multispecies coalescent, it is assumed that all genes from a same species are orthologous, neglecting the possibility of duplication, while duplication and loss models assume that sequences mapped to one species are necessarily the product of a duplication (Rasmussen and Kellis 2012).

Recently a maximum likelihood (ML) supertree approach has been proposed, such that the probabilities are associated to errors in the gene trees (M. Steel and Rodrigo 2008; Cotton and Wilkinson 2009). That is, that given the supertree the probability will reflect a penalty against incongruent input trees, regardless of the evolutionary process responsible for the incongruence. The model assumes that a true supertree generates trees with error, such that their probabilities decrease exponentially according to the distance from the original supertree. The incorporation of supertree approaches into a likelihood framework allows us to compare alternative solutions statistically and to be more explicit about our assumptions, incorporating model selection and hypothesis testing. The probability $P(T_i\mid\tau,\lambda)$ of supertree $\tau$ generating tree $T_i$ is given by

$$
%        P(T_i\mid\tau,\lambda) = \frac{e^{-d(T_i,\tau)/\lambda}}{Z(\tau,\lambda)}
        P(T_i\mid\tau,\lambda) \propto e^{-d(T_i,\tau)/\lambda}
$$

where the scale parameter $\lambda$ controls how strongly discordant trees are penalized. The distance $d()$ is usually a discrete tree-to-tree metric that does not take branch lengths into account, and therefore the original ML supertree approach is based on tree topologies only. Furthermore, it was always assumed that each input tree contained a subset of the species, as is commonplace for supertree methods. The particular distance employed will depend on the problem at hand, as e.g. the Robinson-Foulds distance will lead to a supertree in agreement with the largest number of clusters (or clades) from the input trees (Mukul S. Bansal et al. 2010). Likewise, the subtree prune-regraft (SPR) distance (Beiko and Hamilton 2006) has been used in an exponential model of phylogenetic recombination, where the penalty is against trees in adjacent alignment segments that can only be reconciled through several recombination events, as estimated by the SPR distance between them (Oliveira Martins, Leal, and Kishino 2008; Oliveira Martins and Kishino 2010).

Probabilistic supertree approaches like ML supetrees have a lot of potential, as they can leverage statistical sophistication with computational tractability, being able to deal with several processes of incongruence at once (Cotton and Wilkinson 2009). In this work we further extended the ML supertree approach in several notable ways, creating a hierarchical Bayesian supertree model (Ronquist, Huelsenbeck, and Britton 2004). Firstly, we use a multivariate description of the disagreement between each input tree and the supertree, such that several distances can be used at once. Importantly, among the distances we use jointly are the most parsimonious reconciliation costs under the DL and the ILS models (Guigó, Muchnik, and Smith 1996; Than and Nakhleh 2009). In our model the input trees represent gene families and the supertree is the species phylogeny. This allows us to work with more general data sets than supertree methods usually handle, since we can work with input trees with several leaves representing the same species (even larger, in fact, than the supertree). Furthermore we incorporate this multivariate distribution of tree-supertree distances into a hierarchical Bayesian model, such that the uncertainty of the gene trees as well as in the scale parameters $\lambda$ are taken into account.

Our full model can include the phylogenetic likelihood describing how gene trees are supported by the alignments, but in practice we adopted an importance sampling step such that any Bayesian phylogenetic program can provide the individual distribution of gene trees. (That is why we call our method a 'supertree', in connection to the ML supertree (M. Steel and Rodrigo 2008), despite the full model can in principle work with sequence alignments as well). The resulting hierarchical Bayesian model accounts then for the uncertainty within each gene family to estimate the posterior distribution of all compatible species trees. Unlike the fully probabilistic models described above, our approach is based on the simplifying assumption that summaries based on the most parsimonious reconciliation scenarios are sufficient to explain the disagreement between gene and species trees, regardless of its causes. Yet we do not rely on known gene trees, as with most supertree approaches. Instead we work with marginal gene tree distributions from independent Bayesian phylogenetic analyses of multilocus alignments. We implemented this model in a program called *guenomu*, that receives as input a set of unrooted gene tree distributions, one per gene family, as well as a list of species names, and returns the posterior distribution of rooted species trees as well as a posterior distribution of gene trees for each gene family.

In this manuscript we evaluate our method under simulated scenarios in which gene trees evolve inside species tree with gene duplications and losses as well as incomplete lineage sorting, and assuming gene tree inference error. Usually, these different sources of gene tree disagreement are considered separatedly, although unrecognized processes can have a deleterious effect on the species tree estimation (Rasmussen and Kellis 2012), since even a single ancestral duplication can mimic the effect of several deep coalescences. We circumvent this by favoring species trees that are compatible with different phenomena at once, as dictated by our multivariate distance distribution. Each gene family, on its turn, is allowed to have different contributions in each regard, since we do not expect the contribution of duplications and losses, for instance, to be uniform along the genome. Our simulations suggest that our Bayesian hierarchical model is able to reduce gene tree uncertainty and results in reasonable species trees, considering the complexity of the simulation scenario. Unexpectedly, we also found that GTP methods can be quite accurate as well, while coalescent summary statistics methods performed badly under the scenarios considered.
## Materials and Methods
### Definitions
We define *gene family* $i$ as the set of homologous sequences that compose an alignment $D_i$, which can be comprised of paralogs and orthologs belonging to the same or distinct individuals from a diverged population, that we will call a *species*. The only information needed is a mapping between each member of this gene family and the species which it belongs to (species are thus defined as the entities whose phylogenies we are interested in). Importantly, we assume that we do not know beforehand which members of a gene family are related through a duplication and which diverged trhough coalescences. That is, we map the leaves from the gene tree directly to the species tree. We will then call *gene tree* any phylogenetic tree connecting all members of a gene family, whose corresponding random variable for gene family $i$ will be represented by $G_i$. We are particularly interested in the posterior distribution of gene trees as e.g. estimated by Bayesian phylogenetic inference programs like MrBayes or PhyloBayes (Ronquist et al. 2012; Lartillot and Philippe 2004). These will be our *input gene tree distributions*, one per gene family, and for each input distribution of gene trees we will have an associated *posterior distribution*, in a manner analogous to the BUCKy model (Ané et al. 2007). Although we will mostly use the term *tree*, the input and output distributions contain only the tree topologies. That is, in the current implementation of our model we do not use the branch length information, unless otherwise stated. We will resort to the terms *phylogeny* and *topology* when we want to be explicit about the presence or absence, respectively, of branch lengths.
### Bayesian hierarchical model for species trees
We have devised a Bayesian hierarchical model in which the posterior distribution of the species tree $S$ and all other parameters $\Theta$ given a set of $N$ gene family alignments $\bm{D}=(D_1,...,D_N)$ can be described by

$$
\label{eq:postprob}
P(S,{\bm \Theta}\mid\bm{D})\propto
P(\bm{D},\bm{\theta}\mid\bm{G}) P(\bm{G}\mid\bm{\lambda},S) P(\bm{\lambda}\mid\bm{\lambda_0}) P(\bm{\lambda_0}) P(S)
$$

 where $G=(G_1,\hdots, G_N)$ are the gene trees and $\bm{\lambda}=(\lambda_{ij})$ is the vector with distance scale parameters $j$ for gene family $i$, which depend on the hyperhyperparameters $\bm{\lambda_0}$. The term $P(\bm{D},\bm{\theta}\mid \bm{G})$ is the joint probability of alignments $D$ and parameters $\bm{\theta}$ related to the substitution model -- as we will see this calculation is delegated to the user, so our model works equally for DNA, codon and protein alignments. All species trees $S$ are equally probable a priori, leading to the uniform prior $P(S)=|\mathcal{S}|^{-1}$ over all possible rooted trees $S\in\mathcal{S}$. Equation [\[eq:postprob\]](#eq:postprob) can be further decomposed into

$$
\begin{aligned}
  P(S,{\bm \Theta}\mid\bm{D})&\propto&
        P(\bm{\lambda_0})P(S) \times \prod_{i=1}^N P(D_i,{\bm \theta_i}\mid G_i)
        P(G_i\mid\lambda_{i\cdot},S)P(\lambda_{i\cdot}\mid\bm{\lambda_0})
\end{aligned}
$$

 where we see that the hyperhyperparameters $\bm{\lambda_0}$ and the species tree $S$ are shared among all gene families, while other parameters are allowed to vary between families.

(M. Steel and Rodrigo 2008) assumed that the probability $P(G_i\mid S)$ of a gene tree $G_i$ given the species tree $S$ followed an exponential distribution that depended solely on the dissimilarity between $G_i$ and $S$, which can be represented by a distance $d(G_I,S)$ and a scale parameter $\lambda_i$. We generalized this distribution and assume that this probability depends on the disagreement between $G_i$ and $S$ with respect to several distances $d(G_i,S)=(d_1(G_i,S),\hdots,d_J(G_i,S))$ where each measure $d_j$ $(j=1,...,J)$ represents a distinct distance. Therefore the distribution can be written as

$$
P(G_i\mid S) =
        \frac{e^{-\sum_{j=1}^J d_j(G_i,S)/(m_{ij}\lambda_{ij})}}{\sum_{g\in\mathcal{G}_i} e^{-\sum_{j=1}^J d_j(g,S)/(m_{ij}\lambda_{ij})}} =
        \frac{e^{-\sum_{j=1}^J d_j(G_i,S)/(m_{ij}\lambda_{ij})}}{Z_i(S,\lambda_{i\cdot})}
$$

where $Z_i(S,\lambda_{i\cdot})$ is the *normalization constant* (also called partition function; see below) over the set $G_i$ of all possible topologies for gene family $i$, and the $m_{ij}$ are fixed constants.

Each gene family $i$ has a vector $\lambda_{i\cdot}=(\lambda_{i1},\hdots,\lambda_{iJ})$ where each parameter $\lambda_{ij}$ is associated to a distance $d_j$. We have explored two groups of distances, which we will call *reconciliation distances* and *non-parametric distances* (although they are not in general proper metrics since they do not satisfy the symmetry condition). Our model allows for different parameterizations depending on which distances we want to consider for a specific analysis. The reconciliation distances are based on the most parsimonious reconciliations between the rooted species tree and gene trees (Page 1994; Maddison 1997), which can be calculated using the least common ancestor (LCA) mapping between each node $g$ of the gene tree and a corresponding node $M(g)$ of the species tree (Zmasek and Eddy 2001; Mukul S. Bansal and Eulenstein 2013). This mapping can thus be used in two independent ways: to find the minimum number of duplications and losses needed to explain the disagreement (Guigó, Muchnik, and Smith 1996; Zhang 1997), or to find the minimum number of deep coalescences necessary to make the gene tree be congruent with the species tree (Than and Nakhleh 2009; T. Wu and Zhang 2011).

Under the duplication-loss reconciliation model and given the mapping $M(\cdot)$ between nodes of the gene tree $G_i$ rooted at $r$, denoted $G_i^r$, and the rooted species tree $S$, a gene tree node $g$ will represent a duplication if $M(g)=M(g_1)$ or $M(g)=M(g_2)$, where $g_1$ and $g_2$ are the two children of node $g$. The number of duplications $\mathrm{DUPS}(G_i^r,S)$ will be the number of such nodes, and the number of losses $\mathrm{LOSS}(G_i^r,S)$ is calculated in postorder based on the LCA mapping and the duplication nodes. Likewise, under the reconciliation coalescent model we assume all gene tree nodes are coalescences, and we associate to each node $s$ on the species tree the number of extra lineages, which is the number of branches from the gene tree passing through $s$ minus one. The number of deep coalescences $\mathrm{ILS}(G_i^r,S)$ -- also called incomplete lineage sorting events -- is the sum of extra lineages over all species tree nodes. Since we are using the LCA mapping these numbers will be minimal over all possible reconciliations -- for more details, see (T. Wu and Zhang 2011; Mukul S. Bansal and Eulenstein 2013). Since we do not know the root location of the gene trees, we apply the LCA mapping and calculate the minimum costs over all possible roots $r$ to define the optimal root location. From the duplication-loss model we will then have $d_1(G_i,S)=\mathrm{DUPS}(G_i^{r^*},S)$ and $d_2(G_i,S)=\mathrm{LOSS}(G_i^{r^*},S)$ which will be the smallest amongst all $r$, and for for the coalescent model we will have $d_3(G_i,S)=\mathrm{ILS}(G_i^{r^{**}},S)$, where the optimal root $r^{**}$ might be the same as $r^{*}$ or not.

As for the non-parametric distances, they can be any estimate of disagreement between rooted or unrooted trees, of which the most popular example is the Robinson-Foulds (RF) (or symmetric) distance (Robinson and Foulds 1981). Other examples include the SPR distance (Beiko and Hamilton 2006) and minimum weight matchings between tree branches (Nye, Liò, and Gilks 2006; Bogdanowicz 2008). They are non-parametric in the sense that they do not try to model the biological reason for the disagreement, only the outcome. The RF distance can be calculated only when each leaf from one tree is mapped to at most one tree from the other, restricting therefore the gene families that can be considered. However very recently a generalization of the RF distance called $mulRF$ was introduced that relaxes this constraint, allowing for one of the trees to have several leaves with the same label (Chaudhary, Burleigh, and Fernández-Baca 2013). Such multi-labelled tree, or mul-tree, is here the gene tree. Although the original RF distance is between unrooted topologies, the RF distance can be used to compare a rooted tree to another one, rooted or unrooted (Day 1985; Górecki and Eulenstein 2012). Both the $RF$ and the $mulRF$ distances are calculated as the sum of bipartitions present in one tree but not in the other, but for the $mulRF$ distance the species tree is extended such that all leaves representing a species labelled more than once in the gene tree are replaced by a multifurcation with same number of leaves as in the gene tree (Chaudhary, Burleigh, and Fernández-Baca 2013). Therefore we have that $d_4(G_i,S)=\mathrm{mulRF}(G_i,S)$, noting that this value is the same for any rooting of $G_i$ and S since we use the unrooted version of the trees. Like the reconciliation distances, the non-parametric ones ignore branch lengths, since they measure the disagreement between the topologies only.

Our hierarchical model is completed by specifying the priors for the distance scale parameters. For each distance $j$, the scale parameter $\lambda_{ij}$ of gene family $i$ will follow an exponential hyperprior distribution

$$
P(\lambda_{ij}\mid \bm\lambda_0)=\frac{e^{-\lambda_{ij}/\lambda_{0j}}}{\lambda_{0j}}
$$

 where the scale parameter $\lambda_{0j}$ is the $j$-th element of the (hyperhyperparameter) vector $\bm\lambda_0$, shared among gene families. As mentioned in (M. Steel and Rodrigo 2008) a naive estimation of the $\lambda_{i\cdot}$ scale parameters will lead to very permissive values, and thus we employ an informative hierarchical exponential model for them -- that is, we explicitly inform the model that we prefer small distances overall. The $\lambda_{0j}$ variables themselves come from exponential distributions with a fixed parameter. Since distinct gene families will usually have a different number of members, we must ensure that scale parameters $\lambda_{ij}$ can be comparable across gene families. This is achieved through the scaling parameters $m_{ij}$, which are equivalent to a standardization of all distances $d_j(\cdot,\cdot)$ to the interval $[0,1]$. The choice of a common prior parameter $\lambda_{0j}$ shared across gene families helps to avoid overparameterization. But more importantly, it assumes that the underlying biological process reflected by the distance $j$ -- like the duplication rate, for instance -- acts on the whole genome, while at the same time the $\lambda_{ij}$ parameters for the observed distances are particular to each gene family $i$.
### Sampling from the posterior distribution
Sampling directly from the above hierarchical model is very hard, so we employ a Markov chain Monte Carlo sampling (J. S. Liu 2001) where variables are updated one by one using Generalized Multiple-try Metropolis (GMTM) updates (J. Liu, Liang, and Wong 2000; Pandolfi, Bartolucci, and Friel 2010). The GMTM algorithm for proposing an update of current state $x$ to a new state $y$ can be described as

1.  Given current state $x$, draw $k$  samples $y_1,\hdots,y_k$ independently from the proposal distribution $p(\cdot\mid~x)$

2.  Calculate the weights $\overline{w}_i(y_i;x)=\frac{w_i(y_i;x)}{\sum_{l=1}^k w_l(y_l;x)}$ where $w_i(y_i;x)$ can be any bounded, positive function

3.  Select $y=y_j\in{y_1,\hdots,y_k}$ according to $\overline{w}_i(y_i;x)\quad i=1,\hdots,k$ and set $W_y=\overline{w}_j(y;x)$ (assuming $j$ is the chosen index)

4.  Similarly, draw $x_i^*\sim p(\cdot\mid~y)$ for $i=1,\hdots, j-1, j+1,\hdots, k$, and let $x_j^*=x$ (that is, using same index $j$ as before). Calculate $\overline{w}_i(x_i^*;y)=\frac{w_i(x_i^*;y)}{\sum_{l=1}^k w_l(x_l^*;y)}$ and set $W_x=\overline{w}_j(x;y)$

5.  Accept $y_j$ as new state with probability

$$
\min\left(1,\frac{\pi(y_j)p(x\mid y_j)}{\pi(x)p(y_j\mid x)}\frac{W_x}{W_y}\right)
$$

where $\pi(y)$ is the posterior distribution under $y$. A common choice for the weight $w(y;x)$ is the posterior distribution $\pi(y)$ (terms independent from $y$ will cancel out). This is how the vector $\bm\lambda_0$ is updated.

When updating the gene trees, as mentioned before, we avoid direct calculation of the joint probability of alignment $D_i$ and model $\bm\theta_i$, $P(\bm{D},\bm{\theta}\mid \bm{G})$, assumed to be independent between gene families, and that can be expanded as the product of the phylogenetic likelihood and the prior of the substitution model (Yang and Rannala 1997; Larget and Simon 1999)

$$
P(\bm{D},\bm{\theta}\mid \bm{G}) = \prod_{i=1}^N P(D_i,{\bm \theta_i}\mid G_i) = \prod_{i=1}^N P(D_i\mid{\bm \theta_i}, G_i) P(\bm \theta_i)
$$

 We furthermore assume that for each gene family $i$ we have a posterior distribution of trees $P(G_i\mid D_i, {\bm\theta_i})$ (for example estimated by some state-of-the-art phylogenetic inference program) and use the posterior frequency of tree $\hat{G}$ as a resampling weight $g_i(\hat{G})= P(G_i=\hat{G}\mid D_i, {\bm\theta_i})\propto P(D_i,{\bm \theta_i}\mid\hat{G})$ (Smith and Gelfand 1992). This procedure is similar to the two-stage MCMC approach implemented in BUCKy (Ané et al. 2007), and in practice means that by using the resampling weights as the proposal distribution $p(\cdot\mid\cdot)$ in the GMTM updates the terms $P(D_i,{\bm \theta_i}\mid G_i)$ will cancel out. Under this approach the substitution parameters $\bm\theta_i$ and the branch lengths are nuisance parameters and are neglected from the analysis. At every iteration we propose a new species tree by applying one or several branch swaps on the current one, and we also have a rerooting proposal that does not change the topological information.
### The normalization constant
When updating the species tree $S$ or the distance scale parameters $\lambda_{i\cdot}$, we must also update the normalization constant -- also called the partition function (Murray, Ghahramani, and MacKay 2012; Chung, Perna, and Ane 2013). It is the denominator of our multivariate exponential distribution, and ensures that this distribution sums up to one over all possible gene topologies with the same leafset. It must be taken into account if we want to interpret the penalties as probabilities, and neglecting it can even change the relative score of species trees under some distances (Bryant and Steel 2008). This function changes with the species tree and $\lambda_{i.}$, and computing it is impractical for gene trees with more than a few leaves. With the GMTM algorithm the weight function is not restricted to probability measures from the model. Hence we could postpone the partition function calculation to step (5) of the GMTM algorithm by employing the unscaled distance distribution as the weight function. However the acceptance probability will depend on the ratio between the partition functions $Z_i(S^*,\lambda_{i\cdot}^*)/Z_i(S,\lambda_{i\cdot})$ which do not cancel out.

Our solution is to use an exchange algorithm, where the proposal of a parameter that affects the partition function is always accompanied by the proposal of an auxiliary gene tree, therefore canceling out the partition functions (Murray, Ghahramani, and MacKay 2012; Liang 2010; Atchade, Lartillot, and Robert 2008; Caimo and Friel 2012). This auxiliary gene tree is an augmented variable which is unrelated to our resampled estimates -- it is sampled through a secondary MCMC given the proposed $(S^*,\lambda_{i\cdot}^*)$ and is discarded immediately after calculating the exchange ratio (Caimo and Friel 2010, 2012).

Namely, given the unscaled distance distribution $f_i(G,S,\lambda_{i\cdot})=e^{-\sum_{j=1}^J d_j(G,S)/(m_{ij}\lambda_{ij})}$ then the joint distance probability can be written as

$$
P(\bm{G}\mid \bm{\lambda},S) =
        \frac{\prod_{j=1}^N f_i(G_i,S,\lambda_{i\cdot})}{\prod_{i=1}^N Z_i(S,\lambda_{i\cdot})} =
        \frac{f(\bm{G},S,\bm{\lambda})}{\prod_{i=1}^N Z_i(S,\lambda_{i\cdot})}
$$

 Then the GMTM exchange algorithm proceeds as following, remembering that each state is composed of the species tree $S$ and the vector $\bm{\lambda}$, but in practice we update just $S$ or just a block of $\lambda_{i\cdot}$ for a gene family $i$ at a time.

1.  Given current state $x$, draw $k$  samples $y_1,\hdots,y_k$ independently from the proposal distribution $p(\cdot\mid~x)$

2.  Calculate the weights $\overline{w}_i(y_i;x)=\frac{w_i(y_i;x)}{\sum_{l=1}^k w_l(y_l;x)}$ where $w_i(y_i;x)=f(\bm{G},S^{(i)},\bm{\lambda}^{(i)})$ is the joint unscaled distribution

3.  Select $y=y_j\in{y_1,\hdots,y_k}$ according to $\overline{w}_i(y_i;x)\quad i=1,\hdots,k$ and set $W_y=\overline{w}_j(y;x)$

4.  Draw $x_i^*\sim p(\cdot\mid~y)$ for $i=1,\hdots, j-1, j+1,\hdots, k$, and let $x_j^*=x$ Calculate $\overline{w}_i(x_i^*;y)=\frac{w_i(x_i^*;y)}{\sum_{l=1}^k w_l(x_l^*;y)}$ and set $W_x=\overline{w}_j(x;y)$

5.  Assuming that the proposed state is $y_j=(S^*,\bm{\lambda}^*)$, then draw a set of auxiliary gene trees $\bm{A}=(A_1,\hdots,A_N)$ from $P(\cdot\mid\bm{\lambda}^*,S^*)$. This can be achieved by MCMC starting at the current $\bm{G}$ (notice that within this MCMC the partition functions cancel out.)

6.  Accept $y_j=(S^*,\bm{\lambda}^*)$ as new state with probability

$$
\min\left(1,\frac{P(\bm{\lambda}^*\mid\bm{\lambda_0})p(x\mid y_j)}{P(\bm{\lambda}\mid\bm{\lambda_0})p(y_j\mid x)}
     \frac{f(\bm{G}\mid\bm{\lambda}^*,S^*)f(\bm{A}\mid\bm{\lambda},S)}{f(\bm{A}\mid\bm{\lambda}^*,S^*)f(\bm{G}\mid\bm{\lambda},S)}
     \frac{W_x}{W_y}\right)
$$

Our preliminary results suggested that point estimates from the posterior distribution of species trees are not affected by neglecting the partition function, at least when using multiple distances as we did. Therefore, to reduce the computational workload here we will show results where this normalization constant has been neglected.
### Maximum likelihood estimation of species trees
We furthermore implemented a maximum likelihood version of our model, where instead of sampling from the posterior distribution $P(S,\bm{\Theta}\mid\bm{D})$ we try to find only the optimal values $\hat{S}$ and $\hat{\bm{\Theta}}$ by incrementally decreasing the temperature $kT$ of the distribution $P(S,\bm{\Theta}\mid\bm{D})^{1/kT}$. Therefore we can use the same algorithms as before, since as the temperature $kT$ decreases only updates that improve the posterior distribution are accepted (Rubenthaler, Rydén, and Wiktorsson 2009). The final state of this modified MCMC sampler is then its simulated annealing estimate, which can be interpreted as a ML supertree estimate under an arbitrary set of distances (of which DLI and DLIR have been employed in the present study, as we will see next.)
### Software implementation: *guenomu*
We implemented our Bayesian supertree model into the program *guenomu* (<http://bitbucket.org/leomrtns/guenomu>), where the user can choose which distances will be taken into account by the multivariate exponential distribution. Below we work with two parameterizations: one considering only the reconciliation distances, which we will call the DLI model since it employs the minimum number of Duplications, Losses and Incomplete lineage sortings; and another adding the $mulRF$ distance, called the DLIR model. The program *guenomu* is parallelized at the gene family level, whereby the communication is minimized through the implementation of shared pseudo-random generation streams (Feng et al. 2003; Feng, Cameron, and Buell 2006). It also allows for implementing other distances without compromising the model.
### Assessment of accuracy using computer simulations
In order to assess the performance of our implementation, we carried out a simulation study to compare *guenomu's* posterior species tree distribution with the true species trees and with other similar supertree methods described below. We generated 12348 replicate data sets using *SimPhy* (Mallo, Oliveira Martins, and Posada 2014). Our simulation scenarios takes into account the birth and death of (paralogous) loci inside a species tree as well as the coalescent process describing the gene tree within each locus.

For each replicate data set, a species tree was generated under the Kingman coalescent process using the software Dendropy (Sukumaran and Holder 2010) such that the number of species was between $6$ and $100$ and the length from the root to tips was between $10^2$ and $10^4$ generations. We assumed that the effective population size was constant and that the species tree length in coalescent units (Leaché and Rannala 2011) was between $0.05$ and $5$. Then we used *SimPhy* to simulate the evolution of a number of gene families inside each species tree, generating gene trees with branch lengths in substitution units. Input parameters for each replicate were sampled stochastically from the distributions summarized in the simulation-settings table below. In short, for each gene family the number of individuals per species varied between $1$ and $10$ while the number of gene families ranged from between $3$ and $50$. Furthermore each gene family differed from the species tree through a birth-and-death process of duplications and losses such that the average number of duplications was between $10^{-3}$ and $4$. In the end, each replicate had gene trees with average sizes ranging from $19$ to $1430$ leaves. Both the generated species tree and each gene family tree are ultrametric, with lengths given in generation times. But we used *SimPhy* to further simulate substitution rate heterogeneities at several levels, by applying gamma-distributed rate multipliers at the branches of each of the simulated trees from the three-tree model, namely the species tree, the locus tree and the gene tree (Mallo, Oliveira Martins, and Posada 2014) . As a result the gene family trees will not follow a molecular clock. Although our model for species tree inference do not take into account this branch length information, it is relevant to simulate gene tree uncertainty (i.e., tree inference error), as we will see below.
### Gene tree uncertainty
The simulation procedure just described provide us with the single, true tree underlying each gene family. However, our model explicity considers gene tree uncertainty, relying as input on a sample of the posterior distribution of trees for each gene family. To obtain these we could have simulated sequences along the true gene trees (Yang 2007; Fletcher and Yang 2009), and then infer posterior gene tree distributions using popular Bayesian phylogenetics programs like MrBayes (Ronquist et al. 2012), PhyloBayes (Lartillot and Philippe 2004) or BEAST (Drummond et al. 2012). Unfortunately this approach would be extremely slow, limiting the scope of our simulations, so we devised a way of directly building gene tree distributions for each gene family in an attempt to recreate gene tree uncertainty.

For each true gene tree we generated a population of topologies that may differ from the original by one or more nearest-neighbor interchanges (Waterman and Smith 1978), assuming a probability $p_T$ of applying branch swaps to each tree, and a probability of branch swap per branch inversely proportional to its length, bounded by the maximum swap probability $p_B$. The latter accounts for the observation that shorter branches are harder to reconstruct, and generally are where most of the uncertainty in the phylogenetic estimation is concentrated. Branch swaps implied at random one of the two nearest neighbors. In this way, a proportion $(1-p_T)$ of resulting trees are expected to be identical to the original gene tree, while the others should differ by at least one branch. Values for $p_B$ larger than 1 are allowed, meaning that even larger branches will have maximum uncertainty. For low levels of uncertainty the most frequent topologies will match the true gene family trees, although for the parameter values employed in our simulations almost all data sets contained at least one most frequent gene tree distinct from the true one. The pseudo-code for applying uncertainty to a gene tree is described in the Appendix.

The general workflow of the simulation study is represented in [Figure 1 in the published article](https://doi.org/10.1093/sysbio/syu082). Each replicate data set was composed by a species tree and a set of 3 -- 50 gene families, where phylogenetic uncertainty around each gene family is represented by a tree collection with up to 160 distinct topologies. As a result of the introduction of uncertainty, in 8% of simulations all gene tree topology distributions showed as most frequent tree the original, true tree. In the remaining 92% at least one distribution had a gene tree more frequent than the true tree. From these, wrong trees had on average a cumulative frequency 17% higher than the true tree, although in all cases the true topology could be found amongst the input gene trees, even if at a low frequency. These gene tree frequencies were assumed by *guenomu* to be the resampling weights (posterior frequencies from independent analysis) of each tree.
### Comparison with other species tree approaches
For each simulated data set two Bayesian *guenomu* runs were conducted, under both the DLI and DLIR parameterizations. Each run generates a (posterior) distribution of rooted species trees, and their equivalent posterior distribution of unrooted trees. We used as point estimates of this distribution both the *maximum a posteriori (MAP) tree* and the 50% majority rule *consensus tree* -- remembering that they contain only the topological information and no branch lengths. Furthermore, under the DLIR model we conducted another run of our model to estimate a *maximum likelihood* species tree using simulated annealing, as described above.

For the same simulated data sets, and for comparison purposes, we also employed other supertree-related models to reconstruct the species tree, namely Gene Tree Parsimony (Mukul S. Bansal and Eulenstein 2013) and coalescent distance matrix methods (Liang Liu et al. 2009), also called coalescent summary statistics methods. Because in this study we consider gene trees evolving inside a specie tree under GDL and ILS, we did not consider species tree methods that can not deal with multi-labelled gene trees, like the supermatrix approaches or consensus methods like those implemented in BUCKy (Ané et al. 2007).

Gene Tree Parsimony aims to find species tree that minimizes its total reconciliation cost with a given collection of gene tree topologies. Here we used the software *iGTP* (Chaudhary et al. 2010), which assumes that the disagreement between gene trees and the species tree is either due to duplications and losses, or to incomplete lineage sorting. An independent analysis was conducted under each assumption. *iGTP* allows for each gene tree to have a weight associated with it, and therefore we used the topology frequencies resulting after the introduction of phylogenetic uncertainty (see above). However, for computational feasibility we limited the input gene trees to the most frequent quartile.

We implemented in-house several summary statistics, coalescent species tree reconstruction methods based on gene-wise distance matrices, as described in (Helmkamp, Jewett, and Rosenberg 2012). This class of methods start by creating a matrix of distances between species for each gene phylogeny. Such matrices are then merged into a single one, from which a species tree is built using a clustering algorithm. The four procedures implemented are GLASS (Mossel and Roch 2008), STEAC (L. Liu et al. 2009), SD (Maddison and Knowles 2006), and MAC (Helmkamp, Jewett, and Rosenberg 2012) methods, which differ on how they handle contradicting distances for the same pair of species. This distinction stems from the fact that when there is more than one individual from the same species in a given gene family, its tree will have more than one path for the same pair of species. In such cases the GLASS and SD methods will choose the minimum distance between the pair as the element gene-wise distance matrix, while the MAC and STEAC methods will take the average between all distances for the same pair of species. However when there are several gene families each matrix might have a distinct distance between a given pair of species. In these case, the GLASS and the MAC methods use the minimum value, while the SD and STEAC use the mean distance accross matrices. The species tree is then reconstructed based on this final distance matrix, where SD and STEAC employ an UPGMA, and GLASS and MAC use a single-linkage clustering algorithm (Helmkamp, Jewett, and Rosenberg 2012).

These distance matrix methods are grounded on the divergence times between species, and thus rely heavily on the individual branch lengths. Moreover, these methods assume a single gene tree for each gene family. However we assume that, as with the previous approaches, there is gene tree uncertainty. To make a fair comparison with *guenomu* and iGTP we generated gene tree distributions as before but only used the most frequent gene topology, assuming that all branches had a length of $1.0$. The latter corresponds to using the path length (number of nodes) between leaves as the pairwise distances (M. A. Steel and Penny 1993). In addition, and for comparison purposes, we also used the true gene trees with the true branch lengths. From each of these methods we obtain a single estimate of the species tree.

We also did a preliminary analysis with the mulRF supertree model (Chaudhary, Burleigh, and Fernández-Baca 2013) as implemented in-house, but pilot experiments showed very bad performance and we did not explored it further. Due to their computational complexity, we did not evaluate sequence-based approaches like \*BEAST (Heled and Drummond 2010) as each simulated gene family had hundreds of taxa, on average. Indeed, it would have been very interesting to compare our model with the recently developed software PhylDog, which implements a probabilistic model to estimate the maximum likelihood species tree given a set gene families (Boussau et al. 2013). Unfortunately this software is very difficult to install and run, and we were not successful in designing a pipeline that could systematically execute PhylDog against our simulated data sets.
### Performance measures
We called *tree accuracy* to the number of splits (bipartitions or branches) from the true unrooted species tree that were successfully recovered in the inferred species tree, divided by the total number of splits on the true unrooted species tree. This measure will go from zero when no true splits are recovered to one when the true species tree topology is completely recovered. We also recorded the proportion of simulations where the true unrooted species tree was perfectly recovered (i.e., the proportion of replicates where tree accuracy was one, or *true tree recovery*.) In addition, from the analyses were we had a sample from the posterior distribution we also calculated the proportion of credible sets (at the 95% level) that contained the true unrooted species, or *true tree coverage*.
### Empirical benchmark
To test our model on real data, we downloaded all maximum likelihood gene family trees from the TreeFam database (Schreiber et al. 2014) and then pruned all data outside the 12 *Drosophila* species subtree. We considered all the gene families with more than 3 species. Within these we identified 4591 single-copy gene families (that is, just orthologs), and from the remaining 2562 gene families with at least one paralog, we arbitrarily selected the 43 largest ones, with number of leaves between 102 and 295. As in the simulations, we added phylogenetic uncertainty around the TreeFam gene families as decribed above, generating a distribution of gene trees for each gene family.
## Results
In the simulation scenarios considered, our Bayesian supertree implementation performed best regardless of its parameterization or tree point estimate, showing median tree accuracies around 0.7 ([Figure 2 in the published article](https://doi.org/10.1093/sysbio/syu082)). Our Bayesian method was however followed very closely by iGTP, which was slightly more accurate when DL costs instead of ILS costs were used. Our ML estimates were almost identical to iGTP. The GLASS method ranked next when given the true gene trees, but performed badly when gene trees were uncertain. SD methods did worse, performing equally well on true and observed gene family trees. By far, the worst methods were STEAC and in particular MAC, even when given the true gene trees.

To better understand the effect of different parameter values used in the simulation we fitted a multiple linear regression model to the tree accuracies obtained with *guenomu*'s consensus tree under the DLIR model ([Figure 3 in the published article](https://doi.org/10.1093/sysbio/syu082)). Several explanatory variables were significantly and positively correlated with tree accuracy (p-value $<10^{-15}$ unless otherwise stated), including the smallest branch length and the sum of species tree branch lengths in coalescent units, the expected number of duplications per gene family on the species tree, the average number of individuals from the same species per gene family, the average number of distinct species per gene family (p-value $<10^{-5}$) and the average gene family size (p-value $<10^{-6}$) . The species tree height (length from root to tips) was also significantly correlated with tree accuracy according to simple linear regression, but this relationship faded away once we included other variables like the smallest branch and the total sum of branch lengths. In any case, the most influential parameter was apparently the number of species, which was negatively correlated with tree accuracy (p-value$<10^{-15}$).

Neither *guenomu* nor *iGTP* were able to recover the true species tree for more than 35 species ([Figure 4 in the published article](https://doi.org/10.1093/sysbio/syu082)). Under guenomu, the true tree could be found within the 95% coverage on data sets with up to 40 species. Nonetheless, even when the true species tree could not be found, *guenomu* found trees consistently more accurate than *iGTP* for large trees ([Figure 5 in the published article](https://doi.org/10.1093/sysbio/syu082)). In these cases the consensus tree from the posterior samples under *guenomu* was in any case the best estimate. For example, for the simulations with more than 60 species, 58% of the species trees reconstructed from *guenomu* found at least 60% of the true splits, against only 51% of the trees inferred by *iGTP* ([Figure 5 in the published article](https://doi.org/10.1093/sysbio/syu082)). 5.7% of trees from guenomu reconstructed more than 80% of the true splits for the same large data sets, against 3.4% when *iGTP* was used.

Finally, we fit another multiple linear regression model over all replicates where the dependent variable was the difference in tree accuracy between *guenomu*'s consensus estimate from the DLIR model and iGTP under the DL cost. This indicated that *guenomu* significantly outperformed *iGTP* for higher levels of gene tree uncertainty as measured by $p_T$ and $p_B$, as well as for lower duplication rates and more and larger gene families (at a p-value lower than $10^{-7}$ each). The improvement of *guenomu* over *iGTP* was only marginally significant (p-value $<0.05$) for larger species trees with shorter branch lengths.

Besides the posterior distribution of species trees, *guenomu* also obtains for each replicate the posterior distribution of gene trees, which correspond to reweighted versions of the input gene tree distributions. Gene tree posterior distributions were always more compact than the input gene tree distributions, with fewer trees having a larger contribution, and this reduction was more pronounced for the DLIR parameterization than for DLI (results not shown).
### TreeFam Drosophila gene families
As for the TreeFam *Drosophila* data, *guenomu* reconstructed the known species tree (Stark et al. 2007) as well as the root location for the data set with few large gene families (43 genes with more than 100 members each). The root location could be found even in the absence of an outgroup. The rooted species tree estimation from this data set was robust to gene tree uncertainty. That is, even after replacing each ML tree by a distribution of topologies as described before, *guenomu* was capable of finding the well-known species tree with posterior probability equal to one.

However, for the large data set composed of single-copy gene families (4591 gene families with up to 12 members each), although the unrooted topology could be found in the absence of gene uncertainty, the rooted location was not correctly inferred. In the presence of gene uncertainty, then even the reconstructed topology missed one or two splits.
## Discussion
### Reconstructing species trees from dynamic gene families with lineage sorting
Most phylogenomic analyses start with the identification of ortholog sequences, after some process of sequence collection and pruning of paralogs (e.g. (Koonin and Wolf 2008; Dunn et al. 2008; Medina, Jones, and Fitzpatrick 2011; Sanderson, McMahon, and Steel 2011; Bininda-Emonds 2011; Williams et al. 2012; Lang, Darling, and Eisen 2013; Romiguier et al. 2013; Salichos and Rokas 2013).) Much less often, alternative species tree analyses focus on gene duplication and loss (Katz et al. 2012), or consider paralogs at all (Holton and Pisani 2010). The common practice of removing data (sequences or whole genes) before analysing the reduced data sets is sometimes needed, in order to conform to model assumptions or due to computational limitations of the methods. This process of data carving is also rooted on the elusive idea of removing uncertainty and improving resolution, which is anathema to the Bayesian paradigm. Here we propose a Bayesian supertree model that does not need to separate homolog from paralogs, that can work with multiple gene families and that is also able to deal with multiple individuals per species. It also minimizes the pre-processing of data to conform to strict standards of resolution, since the model can cope with uncertainty in the phylogenetic inference. Therefore our approach, implemented in the program *guenomu*, is able to work with all the available data at hand.

Our simulations suggest that *guenomu* is able to reconstruct reasonably accurate species trees in the presence of gene duplications, losses and deep coalescences, even though it uses simple distances that do not take into account branch lengths. Although it is not a panacea, it suggests that exploring these apparently naïve models of disagreement may be an attractive alternative to biologically realistic frameworks, specially once we aim for a consilience of methods. We hope that more analyses are conducted under broadly usable frameworks like the one we offer, at least as a starting point before embarking into other ones with additional assumptions.

Both DLI and DLIR parameterizations performed almost identically in recovering the species tree, but the DLIR model, which includes the mulRF distance (Chaudhary, Burleigh, and Fernández-Baca 2013), helped reducing the uncertainty from the posterior gene tree distributions. Since the computational overhead of the DLIR parameterization over DLI is minimal, it seems worthy to use DLIR.

In the simulations the true species tree became more and more difficult to recover as the number of species increased. However, even in these cases the majority of true splits were recovered through the consensus or MAP estimates, suggesting that our model can be applied to relatively large species trees. We notice however that for the largest species trees the true species tree was not covered in the 95% credible interval, a behaviour expected from most phylogenetic methods.

Surprisingly, iGTP also provided good species tree estimates, being only marginally worse than *guenomu* as the number of species increased. We did not expect this behaviour because iGTP cannot handle both DL and ILS costs at once, and nonetheless assuming just a DL cost seemed to work quite well even in the presence of ILS. Indeed, iGTP is faster than *guenomu*, which then makes it a very competitive approach. On the other hand, *guenomu* can replicate completely *iGTP*'s functionality through its ML algorithm, which could furthermore work with other complex distance costs in the future.

The distance matrix methods evaluated here assume that the source of disagreement comes from the coalescent process, therefore it comes at no surprise that they did not perform well under our simulations, with an important contribution from duplications and losses. Duplications and losses, when unnacounted for, will have a deleterious effect on the estimated coalescent process (Rasmussen and Kellis 2012). Although this comparison might seem somehow unfair to these summary statistics methods, it serves also as a warning against using attractive algorithms without considering its assumptions. It might be interesting, however, to explore more extensively these methods and see how they fare against *guenomu* or *iGTP* when there is only ILS.

We also observed that both the GTP and our multivariate distance models benefit from large gene families, with several representatives from each species -- even though we do not distinguish paralogs or populational samples. This testifies against the trend of finding good, hand-crafted sets of orthologs at the expense of the paralogs that must be discarded from the analysis. This pruning may be done after the complete data set was used to create the posterior distribution, in a subsequent analysis, but not before.

As expected, the length in coalescent units of the smallest branch of the species tree, and the overall species tree length had a great influence on species tree accuracy according to our model and also for GTP methods, together with the number of species. Indeed, we expect more incomplete lineage sorting and consequently more disagreement between gene and species trees next to short branches, and also over a short generation span within large effective population sizes. Under the conditions simulated, our model performed better on data sets with more and larger gene families, but representing fewer species. Small gene families over a large number of species proved more challenging to the methods compared.

The root location of the *Drosophila* largest gene families was found even without an outgroup since the presence of paralogs is informative about the optimal location of duplication nodes (Katz et al. 2012). For the single-copy genes data set, on the other hand, the root location was not correctly inferred for the same reason: orthologous data sets lack such phylogenetic information. Furthermore the smaller gene families were more sensitive to the effect of tree uncertainty.
### Advantages and shortcomings of our Bayesian supertree method
Bayesian methods in general give us information about the reliability of their findings through the posterior distribution. The flatness of the posterior distribution of species trees, or the frequency of the most frequent topology can give us a hint about the amount of signal contained in the data. The use of a Bayesian framework not only allows for the test of hypothesis and for model selection, but also facilitates the creation of models of increased complexity, subject always to computational constraints. Thus our framework can be extended in a straightforward manner to handle more complex situations, perhaps incorporating more biological realism. Our distinct parameterizations allow us to compare models that use different sets of distances, such that we can, in principle, compare their underlying assumptions.

Our modelling of \"phylogenomic error\" can be seen as a generalization of optimal Bayes estimators, like the consensus tree or quartet puzzling (Huggins et al. 2011). Such a model has already been suggested as an alternative to the detailed description of the processes causing the tree incongruences, while at the same time maintaining the convenience of supertree methods (M. Steel and Rodrigo 2008; Cotton and Wilkinson 2009). It can help bridge the gap between sequence-based alignment and phylogenomic analyses (Cotton and Wilkinson 2009), and ours is the first implementation of this class of methods.

Originally, parsimonious reconciliation methods under the DL or the ILS costs worked with rooted gene trees and rooted species trees (Guigó, Muchnik, and Smith 1996; Than and Nakhleh 2009). However, since we define the D, L and ILS distances as the minimum cost over all possible gene root locations, our model works with unrooted input gene trees. This is an important advantage, as in general there is no reliable information about the rooting at the gene tree level, and fixing a root location can severely influence the distances.

The current implementation of our model uses only the topological information, ignoring branch lengths, but it can be extended in the future to deal explicitly with branch lengths on the gene and species trees. However, it is possible to use crudely the branch lengths as an indication of confidence in the bipartition, as implemented in our generator of gene tree uncertainty.

Other topological distances can be easily implemented as well, taking advantage of our solution to the normalization constant. An example of a promising distance would be one based on the recently developed DLCpar algorithm (Y.-C. Wu et al. 2013), which can find the most parsimonious reconciliation scenario by considering duplications, losses, and deep coalescences at once. Another option would be based on algorithms for finding the minimum cost reconciliation in the presence of duplications, losses and HGT (Tofigh, Hallett, and Lagergren 2011; Doyon et al. 2011). One drawback is that many of these algorithms do not return a single solution, since under these scenarios the best reconciliations are not unique.

When using resampling weights it is important to have a large sample with all reasonable trees represented for each gene family (Smith and Gelfand 1992), specially since even lower likelihood trees might contribute more to the species tree than the ML estimate (Boussau et al. 2013). Therefore ideally under our importance sampling algorithm we should avoid using point estimates of the gene trees (like the ML gene tree, for instance), although for the TreeFam data set it did not hamper the analysis. Nonetheless this simplification allows us to focus on the distribution of tree distances while using well-established software to estimate the individual gene tree distributions, saving time and effort.
### Acknowledgments
We want to thank Nicolas Lartillot for his suggestions on the solution for the partition function. This research was financially supported by the European Research Council (ERC- 2007-Stg 203161-PHYGENOM). DM was also supported by an FPI fellowship (BES-2010-031014) from the Spanish Government at the University of Vigo.

## Appendix
### Pseudocode for the generation of tree uncertainty
Here are the steps to apply uncertainty to a gene tree $T$:

1.  find maximum branch length $t_{MAX}=MAX_i(t_i)$ over all internal branches $i$ of $T$, where $t_i$ is the branch length;

2.  associate to each internal branch $i$ a probability of swap $p_i=p_B\times (1 - t_i/t_{MAX})$;

3.  with probability $(1-p_T)$, create a copy of $T$ unchanged;

4.  Otherwise, with probability $p_T$, apply uncertainty as follows. For each internal branch $i$ in postorder (that is, closer to leaves first), using arbitrary root location and assuming left and right children are labelled $c1$ and $c2$:

    1.  If $c1$ is internal but $c2$ is not, then

        -   with probability $p_{c1}/3$ attach the left child of $c1$ to $i$;

        -   otherwise with probability $p_{c1}/3$ attach the right child of $c1$ to $i$;

        -   otherwise do nothing;

    2.  If $c2$ is internal but $c1$ is not, then

        -   with probability $p_{c2}/3$ attach the left child of $c2$ to $i$;

        -   otherwise with probability $p_{c2}/3$ attach the right child of $c1$ to $i$;

        -   otherwise do nothing;

    3.  If both $c1$ and $c2$ are internal branches, then do one of the following:

        -   with probability $p_{c1}(1-p_{c2})$ attach the left or right child of $c1$ to $i$, at random;

        -   with probability $p_{c2}(1-p_{c1})$ attach the left or right child of $c2$ to $i$, at random;

        -   with probability $p_{c1}p_{c2}/2$ attach the left or right child of $c1$ to $i$, then attach the left or right child of $c2$ to one of $i$, the left or right children of $c1$;

        -   with probability $p_{c1}p_{c2}/2$ attach the left or right child of $c2$ to $i$, then attach the left or right child of $c1$ to one of $i$, the left or right children of $c2$;

        -   with probability $(1-p_{c1})(1-p_{c2})$, do nothing;

If the gene tree does not have branch lengths, then a common value is assumed for all branches. Furthermore if there are several trees in a file, the above procedure is repeated for each one of them, keeping track of their frequencies.

### Simulation settings

|                                                                                                  | **Description**                           | **Symbol** | **Distribution**                         |
|:------------------------------------------------------------------------------------------------:|:------------------------------------------|:----------:|:-----------------------------------------|
|                                       **Species tree**                                       |                                           |            |                                          |
|                                           *(Dendropy)*                                           | Number of species                         |            | Uniform(6,$10^2$)                        |
|                                                                                                  | Number of generations (total tree height) |            | Uniform($10^2$,$10^4$)                   |
|                                                                                                  | Expected number of duplications           | $E_{dup}$  | Uniform($10^{-3}$, 4)                    |
|                                                                                                  | Number of gene families                   |    $N$     | Uniform(3,50)                            |
|                                                                                                  | Rate heterogeneity multiplier             |   $H_s$    | Gamma(1,1)                               |
|                                     **Locus tree**                                     |                                           |            |                                          |
|                                            *(SimPhy)*                                            | Gene duplication rate ${}^b$              |  $\beta$   | Exponential($E_{dup}/\sigma$)            |
|                                                                                                  | Gene loss rate                            |            | Uniform(0, $0.75\times\beta$)            |
|                                                                                                  | Rate heterogeneity multiplier             |   $H_l$    | Gamma(1,1)                               |
|                                     **Gene tree${}^a$**                                      |                                           |            |                                          |
|                                            *(SimPhy)*                                            | Effective population size                 |   $N_e$    | $2000 \times$ LogNormal(0,0.25)          |
|                                                                                                  | Number of individuals per species ${}^c$  |            | Uniform(1,10)                            |
|                                                                                                  | Substitution rate (per time unit)         |            | $0.001$                                  |
|                                                                                                  | Rate heterogeneity multiplier             |   $H_g$    | Gamma(1,1)                               |
|                                        **Gene tree**                                         |                                           |            |                                          |
|                                         **uncertainty**                                          | Maximum number of generated trees         |            | 160                                      |
|                                            (inhouse)                                             | Tree dispersion term                      |   $L_T$    | Uniform(2,5)                             |
|                                                                                                  | Tree location term                        |   $D_T$    | Uniform(3,6)                             |
|                                                                                                  | Frequency of trees with uncertainty       |   $p_T$    | $1.1 \times$ Beta($L_T\times D_T$,$D_T$) |
|                                                                                                  | Branch location term                      |   $L_B$    | Uniform(1,5)                             |
|                                                                                                  | Branch-wise error probability             |   $p_B$    | $1.5 \times$ Beta($L_B$,1)               |
|                               ${}^a$ terminology used by *SimPhy*                                |                                           |            |                                          |
| ${}^b$ the term $\sigma$ is the total sum of branch lengths (in generations) in the species tree |                                           |            |                                          |
|          ${}^c$ common for all gene families, emulating perfect sampling of individuals          |                                           |            |                                          |

Parameter values used in the simulations. Each simulation replicate was parameterized with values sampled from predefined statistical distributions. Given a species tree, simulated with *Dendropy*, the program *SimPhy* was used to simulate gene family trees. An inhouse program was then used to mimic gene tree uncertainty, transforming each gene tree in a distribution of gene tree topologies in an attempt to emulate the effect of alignment-based phylogenetic inference in practice.

## References
Akerborg, Orjan, Bengt Sennblad, Lars Arvestad, and Jens Lagergren. 2009. "Simultaneous Bayesian gene tree reconstruction and reconciliation analysis." *Proceedings of the National Academy of Sciences of the United States of America* 106 (14): 5714--19. <https://doi.org/10.1073/pnas.0806251106>.

Ané, Cécile, Bret Larget, David A Baum, Stacey D Smith, and Antonis Rokas. 2007. "Bayesian Estimation of Concordance Among Gene Trees." *Molecular Biology and Evolution* 24 (2): 412--26.

Arvestad, L., a.-C. Berglund, J. Lagergren, and B. Sennblad. 2003. "Bayesian gene/species tree reconciliation and orthology analysis using MCMC." *Bioinformatics* 19 (Suppl 1): i7--15. <https://doi.org/10.1093/bioinformatics/btg1000>.

Atchade, Yves, Nicolas Lartillot, and Christian P. Robert. 2008. "Bayesian Computation for Statistical Models with Intractable Normalizing Constants." *Arxiv* 0804.3152 (April). <https://arxiv.org/abs/0804.3152>.

Bansal, Mukul S, J Gordon Burleigh, Oliver Eulenstein, and David Fernández-Baca. 2010. "Robinson-Foulds supertrees." *Algorithms for Molecular Biology* 5 (January): 18. <https://doi.org/10.1186/1748-7188-5-18>.

Bansal, Mukul S., and Oliver Eulenstein. 2013. "Algorithms for Genome-Scale Phylogenetics Using Gene Tree Parsimony." *Computational Biology and Bioinformatics, IEEE/ACM Transactions on* 10 (4): 939--56.

Beiko, Robert G, and Nicholas Hamilton. 2006. "Phylogenetic identification of lateral genetic transfer events." *BMC Evolutionary Biology* 6 (January): 15. <https://doi.org/10.1186/1471-2148-6-15>.

Bininda-Emonds, Olaf R P. 2004. "The evolution of supertrees." *Trends in Ecology & Evolution* 19 (6): 315--22. <https://doi.org/10.1016/j.tree.2004.03.015>.

---------. 2011. "Inferring the Tree of Life: chopping a phylogenomic problem down to size?" *BMC Biology* 9 ([Figure 2 in the published article](https://doi.org/10.1093/sysbio/syu082)): 59. <https://doi.org/10.1186/1741-7007-9-59>.

Bogdanowicz, Damian. 2008. "Comparing Phylogenetic Trees Using a Minimum Weight Perfect Matching." *2008 1st International Conference on Information Technology*, May, 1--4. <https://doi.org/10.1109/INFTECH.2008.4621680>.

Boussau, Bastien, Gergely J Szöllosi, Laurent Duret, Manolo Gouy, Eric Tannier, and Vincent Daubin. 2013. "Genome-Scale Coestimation of Species and Gene Trees." *Genome Res.* 23 (2): 323--30. <https://doi.org/10.1101/gr.141978.112>.

Bryant, David, and Mike Steel. 2008. "Computing the Distribution of a Tree Metric." *Discrete Mathematics* 6 (3): 420--26. <https://doi.org/10.1109/TCBB.2009.32>.

Caimo, Alberto, and Nial Friel. 2010. "Bayesian inference for exponential random graph models." *Arxiv* 1007.5192 (July): 29. <http://arxiv.org/abs/1007.5192>.

---------. 2012. "Bayesian model selection for exponential random graph models." *Arxiv* 1201.2337 (January): 23. <http://arxiv.org/abs/1201.2337>.

Chaudhary, Ruchi, Mukul S Bansal, André Wehe, David Fernández-Baca, and Oliver Eulenstein. 2010. "iGTP: a software package for large-scale gene tree parsimony analysis." *BMC Bioinformatics* 11 (1): 574. <https://doi.org/10.1186/1471-2105-11-574>.

Chaudhary, Ruchi, J. Gordon Burleigh, and David Fernández-Baca. 2013. "Inferring species trees from incongruent multi-copy gene trees using the Robinson-Foulds distance." *Algorithms for Molecular Biology* 8 (1): 28. <https://doi.org/10.1186/1748-7188-8-28>.

Chung, Yujin, Nicole T. Perna, and Cecile Ane. 2013. "Computing the Joint Distribution of Tree Shape and Tree Distance for Gene Tree Inference and Recombination Detection." *IEEE/ACM Transactions on Computational Biology and Bioinformatics*. <https://doi.org/10.1109/TCBB.2013.109>.

Cotton, James A, and M. Wilkinson. 2009. "Supertrees join the mainstream of phylogenetics." *Trends in Ecology & Evolution* 24 (1): 1--3. <http://linkinghub.elsevier.com/retrieve/pii/S0169534708003261>.

Day, WHE. 1985. "Optimal algorithms for comparing trees with labeled leaves." *Journal of Classification* 2: 7--28. <http://www.springerlink.com/index/q5906x80g44p44k8.pdf>.

Degiorgio, Michael, and James H Degnan. 2010. "Fast and consistent estimation of species trees using supermatrix rooted triples." *Molecular Biology and Evolution* 27 (3): 552--69. <http://www.ncbi.nlm.nih.gov/pubmed/19833741>.

Doyon, Jean-Philippe, Celine Scornavacca, K Yu Gorbunov, Gergely J Szöllősi, Vincent Ranwez, and Vincent Berry. 2011. "An efficient algorithm for gene/species trees parsimonious reconciliation with losses, duplications and transfers." *Comparative Genomics*, no. October: 93--108. <http://www.springerlink.com/index/D8H78633P3N31268.pdf>.

Drummond, Alexei J., Marc A. Suchard, Dong Xie, and Andrew Rambaut. 2012. "Bayesian phylogenetics with BEAUti and the BEAST 1.7." *Molecular Biology and Evolution* 29 (8): 1969--73. <https://doi.org/10.1093/molbev/mss075>.

Dunn, Casey W, Andreas Hejnol, D. Q. Matus, K. Pang, W. E. Browne, S. A. Smith, E. Seaver, et al. 2008. "Broad phylogenomic sampling improves resolution of the animal tree of life." *Nature* 452 (7188): 745--49. <http://www.nature.com/nature/journal/vaop/ncurrent/abs/nature06614.html>.

Edwards, Scott V, Liang Liu, and Dennis K Pearl. 2007. "High-resolution species trees without concatenation." *Proceedings of the National Academy of Sciences of the United States of America* 104 (14): 5936--41. <https://doi.org/10.1073/pnas.0607004104>.

Feng, Xizhou, Duncan A. Buell, John R. Rose, and Peter J. Waddell. 2003. "Parallel Algorithms for Bayesian Phylogenetic Inference." *Journal of Parallel and Distributed Computing* 63 (7-8): 707--18. <https://doi.org/10.1016/S0743-7315(03)00079-0>.

Feng, Xizhou, Kirk W. Cameron, and Duncan A. Buell. 2006. "PBPI: A High Performance Implementation of Bayesian Phylogenetic Inference." *Proceedings of the 2006 ACM/IEEE Conference on Supercomputing*. <https://doi.org/10.1145/1188455.1188535>.

Fletcher, William, and Ziheng Yang. 2009. "INDELible: a flexible simulator of biological sequence evolution." *Molecular Biology and Evolution* 26 (8): 1879--88. <https://doi.org/10.1093/molbev/msp098>.

Górecki, Pawe, and Oliver Eulenstein. 2012. "A Robinson-Foulds Measure to Compare Unrooted Trees with Rooted Trees." *Bioinformatics Research and Applications* 7292: 115--26. [https://doi.org/10.1007/978-3-642-30191-9\\\_12](https://doi.org/10.1007/978-3-642-30191-9\_12).

Guigó, R, I Muchnik, and T F Smith. 1996. "Reconstruction of ancient molecular phylogeny." *Molecular Phylogenetics and Evolution* 6 (2): 189--213. <https://doi.org/10.1006/mpev.1996.0071>.

Heled, Joseph, and Alexei J Drummond. 2010. "Bayesian inference of species trees from multilocus data." *Molecular Biology and Evolution* 27 (3): 570--80. <https://doi.org/10.1093/molbev/msp274>.

Helmkamp, Laura J, Ethan M Jewett, and Noah A Rosenberg. 2012. "Improvements to a Class of Distance Matrix Methods for Inferring Species Trees from Gene Trees." *J. Comput. Biol.* 19 (6): 632--49. <https://doi.org/10.1089/cmb.2012.0042>.

Holton, Thérèse a, and Davide Pisani. 2010. "Deep genomic-scale analyses of the metazoa reject Coelomata: evidence from single- and multigene families analyzed under a supertree and supermatrix paradigm." *Genome Biology and Evolution* 2 (January): 310--24. <https://doi.org/10.1093/gbe/evq016>.

Huggins, P M, W Li, D Haws, T Friedrich, J Liu, and R Yoshida. 2011. "Bayes estimators for phylogenetic reconstruction." *Systematic Biology* 60 (4): 528--40. <https://doi.org/10.1093/sysbio/syr021>.

Katz, L. A., J. R. Grant, L. W. Parfrey, and J. G. Burleigh. 2012. "Turning the Crown Upside down: Gene Tree Parsimony Roots the Eukaryotic Tree of Life." *Systematic Biology* 61 (4): 653--60. <https://doi.org/10.1093/sysbio/sys026>.

Koonin, Eugene V, and Yuri I Wolf. 2008. "Genomics of bacteria and archaea: the emerging dynamic view of the prokaryotic world." *Nucleic Acids Research* 36 (21): 6688--6719. <https://doi.org/10.1093/nar/gkn668>.

Kubatko, Laura Salter, Bryan C Carstens, and L Lacey Knowles. 2009. "STEM: species tree estimation using maximum likelihood for gene trees under coalescence." *Bioinformatics* 25 (7): 971--73. <https://doi.org/10.1093/bioinformatics/btp079>.

Lang, Jenna Morgan, Aaron E. Darling, and Jonathan A. Eisen. 2013. "Phylogeny of Bacterial and Archaeal Genomes Using Conserved Genes: Supertrees and Supermatrices." Edited by Paul J.Editor Planet. *PLoS ONE* 8 (4): e62510. <https://doi.org/10.1371/journal.pone.0062510>.

Lanier, Hayley C., Huateng Huang, and L. Lacey Knowles. 2014. "How Low Can You Go? The Effects of Mutation Rate on the Accuracy of Species-Tree Estimation." *Molecular Phylogenetics and Evolution* 70 (September): 112--19. <https://doi.org/10.1016/j.ympev.2013.09.006>.

Larget, B, and DL Simon. 1999. "Markov chain Monte Carlo algorithms for the Bayesian analysis of phylogenetic trees." *Molecular Biology and Evolution* 16 (6): 750--59. <http://www.renyi.hu/conferences/bp2008/Slides/larget1-budapest.pdf>.

Lartillot, Nicolas, and Hervé Philippe. 2004. "A Bayesian mixture model for across-site heterogeneities in the amino-acid replacement process." *Molecular Biology and Evolution* 21 (6): 1095--1109. <https://doi.org/10.1093/molbev/msh112>.

Leaché, Adam D, and Bruce Rannala. 2011. "The accuracy of species tree estimation under simulation: a comparison of methods." *Systematic Biology* 60 (2): 126--37. <https://doi.org/10.1093/sysbio/syq073>.

Liang, Faming. 2010. "A double Metropolis-Hastings sampler for spatial models with intractable normalizing constants." *Journal of Statistical Computation and Simulation* 80 (9): 1007--22. <https://doi.org/10.1080/00949650902882162>.

Liu, JS, F Liang, and WH Wong. 2000. "The multiple-try method and local optimization in Metropolis sampling." *Journal of the American Statistical Association* 95 (449): 121. <http://www.tandfonline.com/doi/abs/10.1080/01621459.2000.10473908>.

Liu, Jun S. 2001. *Monte Carlo Strategies in Scientific Computing*. Springer.

Liu, Liang, and Lili Yu. 2010. "Phybase: an R package for species tree analysis." *Bioinformatics* 26 (7): 962--63. <https://doi.org/10.1093/bioinformatics/btq062>.

Liu, Liang, Lili Yu, Laura Salter Kubatko, Dennis K Pearl, and Scott V Edwards. 2009. "Coalescent methods for estimating phylogenetic trees." *Molecular Phylogenetics and Evolution* 53 (1): 320--28. <https://doi.org/10.1016/j.ympev.2009.05.033>.

Liu, L., Lili Yu, Dennis K Pearl, and Scott V Edwards. 2009. "Estimating Species Phylogenies Using Coalescence Times among Sequences." *Systematic Biology* 58 (5): 468--77. <https://doi.org/10.1093/sysbio/syp031>.

Maddison, Wayne P. 1997. "Gene Trees in Species Trees." *Systematic Biology* 46 (3): 523--36. <https://doi.org/10.1093/sysbio/46.3.523>.

Maddison, Wayne P, and L Lacey Knowles. 2006. "Inferring phylogeny despite incomplete lineage sorting." *Systematic Biology* 55 (1): 21--30. <https://doi.org/10.1080/10635150500354928>.

Mallo, Diego, Leonardo de Oliveira Martins, and David Posada. 2014. "SimPhy: Comprehensive Simulation of Gene, Locus and Species Trees at the Genome-Wide Level." *Unpublished*. <https://bitbucket.org/dmallo/stidsim>.

Medina, Edgar M, Gary W Jones, and David A Fitzpatrick. 2011. "Reconstructing the fungal tree of life using phylogenomics and a preliminary investigation of the distribution of yeast prion-like proteins in the fungal kingdom." *Journal of Molecular Evolution* 73 (3-4): 116--33. <https://doi.org/10.1007/s00239-011-9461-4>.

Meng, Chen, and Laura Salter Kubatko. 2009. "Detecting hybrid speciation in the presence of incomplete lineage sorting using gene tree incongruence: a model." *Theoretical Population Biology* 75 (1): 35--45. <https://doi.org/10.1016/j.tpb.2008.10.004>.

Mossel, Elchanan, and Sebastien Roch. 2008. "Incomplete lineage sorting: consistent phylogeny estimation from multiple loci." *IEEE/ACM Transactions on Computational Biology and Bioinformatics / IEEE, ACM* 7 (1): 166--71. <https://doi.org/10.1109/TCBB.2008.66>.

Murray, Iain, Zoubin Ghahramani, and David MacKay. 2012. "MCMC for doubly-intractable distributions." *Arxiv* 1206.6848 (June). <http://arxiv.org/abs/1206.6848>.

Nye, Tom M W, Pietro Liò, and Walter R Gilks. 2006. "A novel algorithm and web-based tool for comparing two alternative phylogenetic trees." *Bioinformatics* 22 (1): 117--19. <https://doi.org/10.1093/bioinformatics/bti720>.

Oliveira Martins, Leonardo de, and Hirohisa Kishino. 2010. "Distribution of Distances Between Topologies and Its Effect on Detection of Phylogenetic Recombination." *Annals of the Institute of Statistical Mathematics* 62 (1): 145--59. <https://doi.org/10.1007/s10463-009-0259-8>.

Oliveira Martins, Leonardo de, Elcio Leal, and Hirohisa Kishino. 2008. "Phylogenetic Detection of Recombination with a Bayesian Prior on the Distance Between Trees." *PLoS ONE* 3 (7): e2651. <https://doi.org/10.1371/journal.pone.0002651>.

Page, RDM. 1994. "Maps Between Trees and Cladistic Analysis of Historical Associations Among Genes, Organisms and Areas." *Systematic Biology* 43 (1): 58--77. <http://sysbio.oxfordjournals.org/content/43/1/58.short>.

Pandolfi, Silvia, Francesco Bartolucci, and PN Friel. 2010. "A generalization of the Multiple-try Metropolis algorithm for Bayesian estimation and model selection." *AISTATS* 9: 581--88. <http://eprints.pascal-network.org/archive/00007723/>.

Ragan, M. A. 1992. "Phylogenetic Inference Based on Matrix Representation of Trees." *Mol Phylogenet Evol* 1 (1): 53--58. <http://www.ncbi.nlm.nih.gov/pubmed/1342924>.

Rannala, Bruce, and Ziheng Yang. 2003. "[Bayes Estimation of Species Divergence Times and Ancestral Population Sizes Using DNA Sequences from Multiple Loci.](https://www.ncbi.nlm.nih.gov/pubmed/12930768)" *Genetics* 164 (4): 1645--56.

Rasmussen, Matthew D, and Manolis Kellis. 2011. "A Bayesian approach for fast and accurate gene tree reconstruction." *Molecular Biology and Evolution* 28 (1): 273--90. <https://doi.org/10.1093/molbev/msq189>.

---------. 2012. "Unified modeling of gene duplication, loss, and coalescence using a locus tree." *Genome Research* 22 (4): 755--65. <https://doi.org/10.1101/gr.123901.111>.

Robinson, DF, and LR Foulds. 1981. "Comparison of phylogenetic trees." *Mathematical Biosciences* 53: 131--47. <http://www.sciencedirect.com/science/article/pii/0025556481900432>.

Romiguier, Jonathan, Vincent Ranwez, Frédéric Delsuc, Nicolas Galtier, and Emmanuel J P Douzery. 2013. "Less Is More in Mammalian Phylogenomics: AT-Rich Genes Minimize Tree Conflicts and Unravel the Root of Placental Mammals." *Mol. Biol. Evol.* 30 (9): 2134--44. <https://doi.org/10.1093/molbev/mst116>.

Ronquist, Fredrik, John P. Huelsenbeck, and Tom Britton. 2004. "Bayesian Supertrees." In *Phylogenetic Supertrees*, edited by Olaf R. P. Bininda-Emonds, 193--224. Springer.

Ronquist, Fredrik, Maxim Teslenko, Paul van der Mark, Daniel L Ayres, Aaron Darling, Sebastian Höhna, Bret Larget, Liang Liu, Marc A Suchard, and John P Huelsenbeck. 2012. "MrBayes 3.2: Efficient Bayesian Phylogenetic Inference and Model Choice Across a Large Model Space." *Systematic Biology*, February, 1--14. <https://doi.org/10.1093/sysbio/sys029>.

Rubenthaler, Sylvain, Tobias Rydén, and Magnus Wiktorsson. 2009. "Fast Simulated Annealing in Rd with an Application to Maximum Likelihood Estimation in State-Space Models." *Stochastic Processes and Their Applications* 119 (6): 1912--31. <https://doi.org/10.1016/j.spa.2008.09.007>.

Salichos, Leonidas, and Antonis Rokas. 2013. "Inferring Ancient Divergences Requires Genes with Strong Phylogenetic Signals." *Nature* 497 (7449): 327--31. <https://doi.org/10.1038/nature12130>.

Sanderson, Michael J, Michelle M McMahon, and Mike Steel. 2011. "Terraces in phylogenetic tree space." *Science (New York, N.Y.)* 333 (6041): 448--50. <https://doi.org/10.1126/science.1206357>.

Schreiber, F., M. Patricio, M. Muffato, M. Pignatelli, and A. Bateman. 2014. "TreeFam V9: A New Website, More Species and Orthology-on-the-Fly." *Nucleic Acids Research* 42 (D1): D922--25. <https://doi.org/10.1093/nar/gkt1055>.

Sjostrand, J., B. Sennblad, L. Arvestad, and J. Lagergren. 2012. "DLRS: gene tree evolution in light of a species tree." *Bioinformatics* 28 (22): 2994--95. <https://doi.org/10.1093/bioinformatics/bts548>.

Smith, Adrian FM, and Alan E Gelfand. 1992. "Bayesian Statistics Without Tears: A Sampling--Resampling Perspective." *The American Statistician* 46 (2): 84--88. <http://www.jstor.org/stable/2684170>.

Stark, Alexander, Michael F. Lin, Pouya Kheradpour, Jakob S. Pedersen, Leopold Parts, Joseph W. Carlson, Madeline A. Crosby, et al. 2007. "Discovery of Functional Elements in 12 Drosophila Genomes Using Evolutionary Signatures." *Nature* 450 (7167): 219--32. <https://doi.org/10.1038/nature06340>.

Steel, M. A., and David Penny. 1993. "Distributions of Tree Comparison Metrics--Some New Results." *Systematic Biology* 42 (2): 126--41. <https://doi.org/10.1093/sysbio/42.2.126>.

Steel, Mike, and Allen Rodrigo. 2008. "Maximum likelihood supertrees." *Systematic Biology* 57 (2): 243--50. <https://doi.org/10.1080/10635150802033014>.

Sukumaran, Jeet, and Mark T Holder. 2010. "DendroPy: A Python Library for Phylogenetic Computing." *Bioinformatics* 26 (12): 1569--71. <https://doi.org/10.1093/bioinformatics/btq228>.

Szollosi, G. J., W. Rosikiewicz, B. Boussau, E. Tannier, and V. Daubin. 2013. "Efficient Exploration of the Space of Reconciled Gene Trees." *Systematic Biology* 62 (6): 901--12. <https://doi.org/10.1093/sysbio/syt054>.

Than, Cuong, and Luay Nakhleh. 2009. "Species tree inference by minimizing deep coalescences." *PLoS Computational Biology* 5 (9): e1000501. <https://doi.org/10.1371/journal.pcbi.1000501>.

Tofigh, Ali, Michael Hallett, and Jens Lagergren. 2011. "Simultaneous identification of duplications and lateral gene transfers." *IEEE/ACM Transactions on Computational Biology and Bioinformatics (TCBB)* 8 (2): 517--35. <https://doi.org/10.1109/TCBB.2010.14>.

Waterman, M. S., and T. F. Smith. 1978. "On the Similarity of Dendrograms." *Journal of Theoretical Biology* 73 (4): 789--800. <https://doi.org/10.1016/0022-5193(78)90137-6>.

Whidden, Chris, Norbert Zeh, and Robert G. Beiko. 2013. "Supertrees Based on the Subtree Prune-and-Regraft Distance." *PeerJ PrePrints* 1 (May): 18. <https://doi.org/10.7287/peerj.preprints.18v1>.

Williams, Tom A, Peter G Foster, Tom M W Nye, Cymon J Cox, and T Martin Embley. 2012. "A congruent phylogenomic signal places eukaryotes within the Archaea." *Proceedings of the The Royal Society / Biological Sciences* 279 (1749): 4870--79. <https://doi.org/10.1098/rspb.2012.1795>.

Wu, Martin, Sourav Chatterji, and Jonathan a. Eisen. 2012. "Accounting For Alignment Uncertainty in Phylogenomics." Edited by Marco Salemi. *PLoS ONE* 7 (1): e30288. <https://doi.org/10.1371/journal.pone.0030288>.

Wu, Taoyang, and Louxin Zhang. 2011. "Structural properties of the reconciliation space and their applications in enumerating nearly-optimal reconciliations between a gene tree and a species tree." *BMC Bioinformatics* 12 Suppl 9 (Suppl 9): S7. <https://doi.org/10.1186/1471-2105-12-S9-S7>.

Wu, Y.-C., M. D. Rasmussen, M. S. Bansal, and M. Kellis. 2013. "Most Parsimonious Reconciliation in the Presence of Gene Duplication, Loss, and Deep Coalescence Using Labeled Coalescent Trees." *Genome Research*, December. <https://doi.org/10.1101/gr.161968.113>.

Yang, Ziheng. 2007. "PAML 4: phylogenetic analysis by maximum likelihood." *Molecular Biology and Evolution* 24 (8): 1586--91. <https://doi.org/10.1093/molbev/msm088>.

Yang, Ziheng, and Bruce Rannala. 1997. "Bayesian phylogenetic inference using DNA sequences: a Markov Chain Monte Carlo Method." *Molecular Biology and Evolution* 14 (7): 717--24. <http://www.ncbi.nlm.nih.gov/pubmed/9214744>.

Yu, Yun, Cuong Than, James H Degnan, and Luay Nakhleh. 2011. "Coalescent histories on phylogenetic networks and detection of hybridization despite incomplete lineage sorting." *Systematic Biology* 60 (2): 138--49. <https://doi.org/10.1093/sysbio/syq084>.

Zhang, L. 1997. "On a Mirkin-Muchnik-Smith conjecture for comparing molecular phylogenies." *Journal of Computational Biology* 4 (2): 177--87. <http://online.liebertpub.com/doi/abs/10.1089/cmb.1997.4.177>.

Zmasek, CM, and SR Eddy. 2001. "A simple algorithm to infer gene duplication and speciation events on a gene tree." *Bioinformatics* 17 (9): 821--28. <http://bioinformatics.oxfordjournals.org/content/17/9/821.short>.
