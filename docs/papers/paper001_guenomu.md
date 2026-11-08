paper https://academic.oup.com/sysbio/article/65/3/397/2468904  https://doi.org/10.1093/sysbio/syu082 _Systematic Biology_, Volume 65, Issue 3, May 2016, Pages 397–416

# A Bayesian Supertree Model for Genome-Wide Species Tree Reconstruction
Reading guide: This is a source-based text extraction, not a summary. Figure images are omitted; captions, mathematics, tables, citations and textual cross-references are retained. The readable article follows. The final “Exact source archive” preserves every supplied LaTeX source and the complete BibTeX database verbatim, including comments, disabled passages, macros and metadata. Consult that archive for exact notation or material excluded from the typeset article; revision notes are not published findings. Figure content available only in image files is outside this text extraction. No external files are needed to read the retained text.
Source document: `paper001_guenomu/20140204.manuscript.tex`. Source files take precedence over any older compiled PDFs in the directory.
or Genome-Wide Species Tree Reconstruction**

[Leonardo de Oliveira Martins$^1$, Diego Mallo$^2$, and David Posada$^3$]{.smallcaps}\
*Department of Biochemistry, Genetics and Immunology, University of Vigo, Vigo, 36310, Spain*\
:::

**Corresponding author:** Leonardo de Oliveira Martins, Department of Biochemistry, Genetics and Immunology, University of Vigo, Vigo, 36310, Spain\
**E-mail:** $^1$leomrtns@uvigo.es, $^2$ dmallo@uvigo.es, $^3$ dposada@uvigo.es\
## Abstract
Current phylogenomic data sets highlight the need for species tree methods able to deal with several sources of gene tree/species tree incongruence. At the same time, we need to make most use of all available data. Most species tree methods deal with single processes of phylogenetic discordance, namely, gene duplication and loss, incomplete lineage sorting or horizontal gene transfer. In this manuscript we address for the first time the problem of species tree inference from multilocus, genome-wide data sets in the presence of gene duplication and loss and incomplete lineage sorting, therefore without the need to identify orthologs or to use a single individual per species. We do this by extending the idea of Maximum Likelihood supertrees to a hierarchical Bayesian model where several sources of gene tree/species tree disagreement can be accounted for in a modular manner. We implemented this model in a computer program called *guenomu* whose input are with posterior distributions of unrooted gene tree topologies for multiple gene families, and whose output is the posterior distribution of rooted species tree topologies. We conducted extensive simulations under complex phylogenomic models in order to evaluate the performance ouf our approach in comparison with other species tree approaches able to deal with multilabeled trees. Our method ranked best under simulated data sets, in spite of ignoring branch lengths, and performed well on empirical data sets as well. Our results show in addition that under complex simulation scenarios , gene tree parsimony is also a competitive approach once we consider its speed in contrast to more sophisticated models.

(Keywords: supertree, reconciliation, hierarchical Bayesian model, tree distance, phylogenomics)\

The evolutionary history of a gene family is not an exact representation of the evolution of the species embedding it, since we must take into account the effects of gene duplication, population factors like incomplete lineage sorting as well as the possibility of horizontal gene transfer (Maddison 1997). In recent years the rapid accumulation of phylogenomic data has highlighted the need for specific methods to infer species phylogenies. Methods of species tree inference based on collection of gene alignments/trees can be broadly classified into supermatrix, supertree and model-based approaches. The supermatrix strategy consists in concatenating different genes into a single large alignment, assuming that minor conflicting phylogenetic signals will cancel out and a common phylogenetic trend -- the species history -- will emerge. Departures from this assumption must be handled by subjectively removing columns or whole regions (M. Wu, Chatterji, and Eisen 2012; Lanier, Huang, and Knowles 2014). This method relies on the correct inference of the orthologous sets, such that the sequences from a single species can be concatenated into a single 'supergene' (Degiorgio and Degnan 2010).

Alternatively, phylogenetic inference can be conducted independently for each gene alignment, and then in a second step the species tree is inferred from the resulting gene trees (Cotton and Wilkinson 2009). In this case the alignments are used indirectly, through summary statistics like maximum likelihood trees, to build a 'supertree'. Under the supertree umbrella, some methods try to find the species tree that minimizes its overall disagreement with the gene tree collection without explicitly taking into account the biological phenomena behind the incongruence, like the RF supertree (Mukul S. Bansal et al. 2010; Chaudhary, Burleigh, and Fernández-Baca 2013) and the MRP approaches (Ragan 1992; Bininda-Emonds 2004). Other supertree approaches aim for the most parsimonious scenario by assuming that the disagreement between the species tree and each of the gene trees is due to HGT (Whidden, Zeh, and Beiko 2013), duplications and losses or deep coalescences (Mukul S. Bansal and Eulenstein 2013). The problem of finding the species tree that minimizes the reconciliation cost is called Gene Tree Parsimony (GTP) (Mukul S. Bansal and Eulenstein 2013).

Most supertree approaches neglect the gene tree branch lengths, but there is another class of species tree methods closely related to the supertree methods which are based on the multispecies coalescent. These methods try to reconstruct the species tree from a matrix of distances between species, which in turn is built based on distance matrices from the individual gene trees (reviewed in (Liang Liu et al. 2009; Helmkamp, Jewett, and Rosenberg 2012)). These methods are based on the observation that the coalescence time of two genes from distinct species always precedes the speciation event. In general they assume that the gene trees are known without error (Kubatko, Carstens, and Knowles 2009; Liang Liu and Yu 2010). Some of these methods also can take into account hibridization (Meng and Kubatko 2009) or HGT (Yu et al. 2011).

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

 where $G=(G_1,\hdots, G_N)$ are the gene trees and $\bm{\lambda}=(\lambda_{ij})$ is the vector with distance scale parameters $j$ for gene family $i$, which depend on the hyperhyperparameters $\bm{\lambda_0}$. The term $P(\bm{D},\bm{\theta}\mid \bm{G})$ is the joint probability of alignments $D$ and parameters $\bm{\theta}$ related to the substitution model -- as we will see this calculation is delegated to the user, so our model works equally for DNA, codon and protein alignments. All species trees $S$ are equally probable a priori, leading to the uniform prior $P(S)=|\mathcal{S}|^{-1}$ over all possible rooted trees $S\in\mathcal{S}$. Equation [\[eq:postprob\]](#eq:postprob){reference-type="ref" reference="eq:postprob"} can be further decomposed into 

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
### The normalization constant {#sec:normal}
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
We implemented our Bayesian supertree model into the program *guenomu* (<http://bitbucket.org/leomrtns/guenomu>), where the user can choose which distances will be taken into account by the multivariate exponential distribution. Below we work with two parameterizations: one considering only the reconciliation distances, which we will call the DLI model since it employs the minimum number of Duplications, Losses and Incomplete lineage sortings; and another adding the $mulRF$ distance, called the DLIR model. The program *guenomu* is paralellized at the gene family level, whereby the communication is minimized through the implementation of shared pseudo-random generation streams (Feng et al. 2003; Feng, Cameron, and Buell 2006). It also allows for implementing other distances without compromising the model.
### Assesment of accuracy using computer simulations
In order to assess the performance of our implementation, we carried out a simulation study to compare *guenomu's* posterior species tree distribution with the true species trees and with other similar supertree methods described below. We generated 12348 replicate data sets using *SimPhy* (Mallo, Oliveira Martins, and Posada 2014). Our simulation scenarios takes into account the birth and death of (paralogous) loci inside a species tree as well as the coalescent process describing the gene tree within each locus.

For each replicate data set, a species tree was generated under the Kingman coalescent process using the software Dendropy (Sukumaran and Holder 2010) such that the number of species was between $6$ and $100$ and the length from the root to tips was between $10^2$ and $10^4$ generations. We assumed that the effective population size was constant and that the species tree length in coalescent units (Leaché and Rannala 2011) was between $0.05$ and $5$. Then we used *SimPhy* to simulate the evolution of a number of gene families inside each species tree, generating gene trees with branch lengths in substitution units. Input parameters for each replicate were sampled stochastically from the distributions summarized in table [1](#tab:params){reference-type="ref" reference="tab:params"}. In short, for each gene family the number of individuals per species varied between $1$ and $10$ while the number of gene families ranged from between $3$ and $50$. Furthermore each gene family differed from the species tree through a birth-and-death process of duplications and losses such that the average number of duplications was between $10^{-3}$ and $4$. In the end, each replicate had gene trees with average sizes ranging from $19$ to $1430$ leaves. Both the generated species tree and each gene family tree are ultrametric, with lenghts given in generation times. But we used *SimPhy* to further simulate substitution rate heterogeneities at several levels, by applying gamma-distributed rate multipliers at the branches of each of the simulated trees from the three-tree model, namely the species tree, the locus tree and the gene tree (Mallo, Oliveira Martins, and Posada 2014) . As a result the gene family trees will not follow a molecular clock. Although our model for species tree inference do not take into account this branch length information, it is relevant to simulate gene tree uncertainty (i.e., tree inference error), as we will see below.
### Gene tree uncertainty
The simulation procedure just described provide us with the single, true tree underlying each gene family. However, our model explicity considers gene tree uncertainty, relying as input on a sample of the posterior distribution of trees for each gene family. To obtain these we could have simulated sequences along the true gene trees (Yang 2007; Fletcher and Yang 2009), and then infer posterior gene tree distributions using popular Bayesian phylogenetics programs like MrBayes (Ronquist et al. 2012), PhyloBayes (Lartillot and Philippe 2004) or BEAST (Drummond et al. 2012). Unfortunately this approach would be extremely slow, limiting the scope of our simulations, so we devised a way of directly building gene tree distributions for each gene family in an attempt to recreate gene tree uncertainty.

For each true gene tree we generated a population of topologies that may differ from the original by one or more nearest-neighbor interchanges (Waterman and Smith 1978), assuming a probability $p_T$ of applying branch swaps to each tree, and a probability of branch swap per branch inversely proportional to its length, bounded by the maximum swap probability $p_B$. The latter accounts for the observation that shorter branches are harder to reconstruct, and generally are where most of the uncertainty in the phylogenetic estimation is concentrated. Branch swaps implied at random one of the two nearest neighbors. In this way, a proportion $(1-p_T)$ of resulting trees are expected to be identical to the original gene tree, while the others should differ by at least one branch. Values for $p_B$ larger than 1 are allowed, meaning that even larger branches will have maximum uncertainty. For low levels of uncertainty the most frequent topologies will match the true gene family trees, although for the parameter values employed in our simulations almost all data sets contained at least one most frequent gene tree distinct from the true one. The pseudo-code for applying uncertainty to a gene tree is described in the Appendix.

The general workflow of the simulation study is represented in Figure [1](#fig:simul){reference-type="ref" reference="fig:simul"}. Each replicate data set was composed by a species tree and a set of 3 -- 50 gene families, where phylogenetic uncertainty around each gene family is represented by a tree collection with up to 160 distinct topologies. As a result of the introduction of uncertainty, in 8% of simulations all gene tree topology distributions showed as most frequent tree the original, true tree. In the remaining 92% at least one distribution had a gene tree more frequent than the true tree. From these, wrong trees had on average a cummulative frequency 17% higher than the true tree, although in all cases the true topology could be found amongst the input gene trees, even if at a low frequency. These gene tree frequencies were assumed by *guenomu* to be the resampling weights (posterior frequencies from independent analysis) of each tree.
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
In the simulation scenarios considered, our Bayesian supertree implementation performed best regardless of its parameterization or tree point estimate, showing median tree accuracies around 0.7 (Fig. [2](#fig:compareall){reference-type="ref" reference="fig:compareall"}). Our Bayesian method was however followed very closely by iGTP, which was slightly more accurate when DL costs instead of ILS costs were used. Our ML estimates were almost identical to iGTP. The GLASS method ranked next when given the true gene trees, but performed badly when gene trees were uncertain. SD methods did worse, performing equally well on true and observed gene family trees. By far, the worst methods were STEAC and in particular MAC, even when given the true gene trees.

To better understand the effect of different parameter values used in the simulation we fitted a multiple linear regression model to the tree accuracies obtained with *guenomu*'s consensus tree under the DLIR model (Fig. [3](#fig:corr){reference-type="ref" reference="fig:corr"}). Several explanatory variables were significantly and positively correlated with tree accuracy (p-value $<10^{-15}$ unless otherwise stated), including the smallest branch length and the sum of species tree branch lengths in coalescent units, the expected number of duplications per gene family on the species tree, the average number of individuals from the same species per gene family, the average number of distinct species per gene family (p-value $<10^{-5}$) and the average gene family size (p-value $<10^{-6}$) . The species tree height (length from root to tips) was also significantly correlated with tree accuracy according to simple linear regression, but this relationship faded away once we included other variables like the smallest branch and the total sum of branch lenghts. In any case, the most influential parameter was apparently the number of species, which was negatively correlated with tree accuracy (p-value$<10^{-15}$).

Neither *guenomu* nor *iGTP* were able to recover the true species tree for more than 35 species (Fig. [4](#fig:true){reference-type="ref" reference="fig:true"}). Under guenomu, the true tree could be found within the 95% coverage on data sets with up to 40 species. Nonetheless, even when the true species tree could not be found, *guenomu* found trees consistently more accurate than *iGTP* for large trees (Fig. [5](#fig:recover){reference-type="ref" reference="fig:recover"}). In these cases the consensus tree from the posterior samples under *guenomu* was in any case the best estimate. For example, for the simulations with more than 60 species, 58% of the species trees reconstructed from *guenomu* found at least 60% of the true splits, against only 51% of the trees inferred by *iGTP* (Fig. [5](#fig:recover){reference-type="ref" reference="fig:recover"}). 5.7% of trees from guenomu reconstructed more than 80% of the true splits for the same large data sets, against 3.4% when *iGTP* was used.

Finally, we fit another multiple linear regression model over all replicates where the dependent variable was the difference in tree accuracy between *guenomu*'s consensus estimate from the DLIR model and iGTP under the DL cost. This indicated that *guenomu* significantly outperformed *iGTP* for higher levels of gene tree uncertainty as measured by $p_T$ and $p_B$, as well as for lower duplication rates and more and larger gene families (at a p-value lower than $10^{-7}$ each). The improvement of *guenomu* over *iGTP* was only marginally significant (p-value $<0.05$) for larger species trees with shorter branch lenghts.

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

::: spacing
1.4
:::
## Appendix {#appendix .unnumbered}
### pseudo-code for the generation of tree uncertainty
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

**Tables**

::: spacing
1.4

::: {#tab:params}
|                                                                                                  | **description**                           | **symbol** | **distribution**                         |
|:------------------------------------------------------------------------------------------------:|:------------------------------------------|:----------:|:-----------------------------------------|
|                                       2-4 **species tree**                                       |                                           |            |                                          |
|                                           *(Dendropy)*                                           | Number of species                         |            | Uniform(6,$10^2$)                        |
|                                                                                                  | Number of generations (total tree height) |            | Uniform($10^2$,$10^4$)                   |
|                                                                                                  | Expected number of duplications           | $E_{dup}$  | Uniform($10^{-3}$, 4)                    |
|                                                                                                  | Number of gene families                   |    $N$     | Uniform(3,50)                            |
|                                                                                                  | Rate heterogeneity multiplier             |   $H_s$    | Gamma(1,1)                               |
|                                     2-4 **locus tree${}^a$**                                     |                                           |            |                                          |
|                                            *(SimPhy)*                                            | Gene duplication rate ${}^b$              |  $\beta$   | Exponential($E_{dup}/\sigma$)            |
|                                                                                                  | Gene loss rate                            |            | Uniform(0, $0.75\times\beta$)            |
|                                                                                                  | Rate heterogeneity multiplier             |   $H_l$    | Gamma(1,1)                               |
|                                     2-4 **gene tree${}^a$**                                      |                                           |            |                                          |
|                                            *(SimPhy)*                                            | Effective population size                 |   $N_e$    | $2000 \times$ LogNormal(0,0.25)          |
|                                                                                                  | Number of individuals per species ${}^c$  |            | Uniform(1,10)                            |
|                                                                                                  | Substitution rate (per time unit)         |            | $0.001$                                  |
|                                                                                                  | Rate heterogeneity multiplier             |   $H_g$    | Gamma(1,1)                               |
|                                        2-4 **gene tree**                                         |                                           |            |                                          |
|                                         **uncertainty**                                          | Maximum number of generated trees         |            | 160                                      |
|                                            (inhouse)                                             | Tree dispersion term                      |   $L_T$    | Uniform(2,5)                             |
|                                                                                                  | Tree location term                        |   $D_T$    | Uniform(3,6)                             |
|                                                                                                  | Frequency of trees with uncertainty       |   $p_T$    | $1.1 \times$ Beta($L_T\times D_T$,$D_T$) |
|                                                                                                  | Branch location term                      |   $L_B$    | Uniform(1,5)                             |
|                                                                                                  | Branch-wise error probability             |   $p_B$    | $1.5 \times$ Beta($L_B$,1)               |
|                               ${}^a$ terminology used by *SimPhy*                                |                                           |            |                                          |
| ${}^b$ the term $\sigma$ is the total sum of branch lengths (in generations) in the species tree |                                           |            |                                          |
|          ${}^c$ common for all gene families, emulating perfect sampling of individuals          |                                           |            |                                          |

: Parameter values used in the simulations. Each simulation replicate was parameterized with values sampled from predefined statistical distributions. Given a species tree, simulated with *Dendropy*, the program *SimPhy* was used to simulate gene family trees. An inhouse program was then used to mimic gene tree uncertainty, transforming each gene tree in a distribution of gene tree topologies in an attempt to emulate the effect of alignment-based phylogenetic inference in practice.
:::
:::

**Figure Captions**

::: {#fig:simul}
**Figure 1 \[fig:simul\] --- image omitted.**

[]{#fig:simul label="fig:simul"} Simulation workflow. On the left a single data set is produced using *Dendropy* and *SimPhy*. The true species tree (rooted, with branch lengths) simulated by *Dendropy* is used as input by *SimPhy* to generate several (rooted) phylogenies, one per gene family. Then, uncertainty is added such that we have a distribution of topologies (unrooted, no branch lengths) per gene family. These collections will be used as input to the inference programs in the right panel. Our software *guenomu* estimates the posterior distribution of species trees (rooted, without branch lengths) and the posterior distributions of gene family topologies, based on all input gene tree distributions. The software *iGTP* also uses the input gene trees after transforming the frequencies into integer values representing weights, and estimates two rooted species trees: one under the Duplication and Loss cost, and another under the Deep Coalescence cost. For the distance-based species tree inference algorithms only one tree per gene family is used, and two alternative choices were attempted: one was to use the true gene families, with branch lengths; and the other was to use the most frequent gene trees (topologies only) after introducing phylogenetic uncertainty. In the later case it was assumed that all branches had the same length.
:::

::: {#fig:compareall}
**Figure 2 \[fig:compareall\] --- image omitted.**

[]{#fig:compareall label="fig:compareall"} Tree accuracy. Each violin plot (kernel density curves plus a boxplot) represents the distribution of tree accuracy values for the different methods evaluated. We ran *guenomu* 3 times for each simulation, to sample from the posterior distribution under the DLI and the DLIR models ('DLI' and 'DLIR'), and to find the ML estimate under the DLIR parameterization ('DLIR-ML'). From each posterior distribution we obtained two point estimates of the species tree: the MAP tree and the consensus tree (labelled 'MAP' and 'cons'). We also obtained two estimates of species trees using iGTP, running the program under the DL cost and under the ILS cost ('iGTP-DL' and 'iGTP-ILS'). We also ran 4 distance matrix approaches (GLASS, SD, STEAC and MAC) using two types of data sets, just the topologies with uncertainty ('topo') or the true simulated gene family trees with branch lengths ('true'). At the top we show the hierarchical clustering of the different methods based on their tree accuracy values for all replicates.
:::

::: {#fig:corr}
**Figure 3 \[fig:corr\] --- image omitted.**

[]{#fig:corr label="fig:corr"} *Guenomu*'s tree accuracy with respect to simulation parameters. Here the species tree was estimated as the consensus tree from the posterior distribution assuming a DLIR parameterization. Based on a multiple linear regression analysis, these are the parameters that most significantly affected tree accuracy: at the top we have, from left to right, the total number of gene families, the total species tree length (in coalescent units), the total number of species on the species tree, and the average gene family size (number of leaves); at the bottom we have the length of the smallest branch in the species tree (in coalescent units), the expected number of duplications per replicate, the average number of individuals from same species per gene family, and the average number of species represented by gene family. Over each panel we show the regression line overlaid in red.
:::

::: {#fig:true}
**Figure 4 \[fig:true\] --- image omitted.**

[]{#fig:true label="fig:true"} Tree accuracy and true tree recovery for *guenomu* and *iGTP* according to the number of species. The lower lines show the true tree recoveries pooled over the number of species, while the higher lines show the average tree accuracy under different conditions. The gray area represents the true tree coverage under *guenomu* DLIR model (ie., the proportion of data sets where the true species tree was within the 95% credible set). Each line is actually a smooth regression line over the individual values. For *iGTP* the species tree was estimated under the DL cost, while for *guenomu* we plot the ML estimate as well as the Bayesian consensus, both under the DLIR parameterization.
:::

::: {#fig:recover}
**Figure 5 \[fig:recover\] --- image omitted.**

[]{#fig:recover label="fig:recover"} Tree accuracies of *iGTP* and *guenomu* for large data sets. On the left we have all data sets where the number of species was between 30 and 60, and on the right we have the pooled replicates with species tree sizes between 61 and 100 leaves. Each bar represents the proportion of data sets with tree accuracies between the ranges described in the middle, such that at the top we have the fraction of simulations where more than 80% of the splits were recovered, while the inferences where less than 20% of the true splits were found are the bottom.
:::
## References
::: {#refs .references .csl-bib-body .hanging-indent}
::: {#ref-Akerborg2009 .csl-entry}
Akerborg, Orjan, Bengt Sennblad, Lars Arvestad, and Jens Lagergren. 2009. "[Simultaneous Bayesian gene tree reconstruction and reconciliation analysis.]{.nocase}" *Proceedings of the National Academy of Sciences of the United States of America* 106 (14): 5714--19. <https://doi.org/10.1073/pnas.0806251106>.
:::

::: {#ref-Ane2007 .csl-entry}
Ané, Cécile, Bret Larget, David A Baum, Stacey D Smith, and Antonis Rokas. 2007. "Bayesian Estimation of Concordance Among Gene Trees." *Molecular Biology and Evolution* 24 (2): 412--26.
:::

::: {#ref-Arvestad2003 .csl-entry}
Arvestad, L., a.-C. Berglund, J. Lagergren, and B. Sennblad. 2003. "[Bayesian gene/species tree reconciliation and orthology analysis using MCMC]{.nocase}." *Bioinformatics* 19 (Suppl 1): i7--15. <https://doi.org/10.1093/bioinformatics/btg1000>.
:::

::: {#ref-Atchade2008 .csl-entry}
Atchade, Yves, Nicolas Lartillot, and Christian P. Robert. 2008. "Bayesian Computation for Statistical Models with Intractable Normalizing Constants." *Arxiv* 0804.3152 (April). <https://arxiv.org/abs/0804.3152>.
:::

::: {#ref-Bansal2010 .csl-entry}
Bansal, Mukul S, J Gordon Burleigh, Oliver Eulenstein, and David Fernández-Baca. 2010. "[Robinson-Foulds supertrees.]{.nocase}" *Algorithms for Molecular Biology* 5 (January): 18. <https://doi.org/10.1186/1748-7188-5-18>.
:::

::: {#ref-Bansal13 .csl-entry}
Bansal, Mukul S., and Oliver Eulenstein. 2013. "Algorithms for Genome-Scale Phylogenetics Using Gene Tree Parsimony." *Computational Biology and Bioinformatics, IEEE/ACM Transactions on* 10 (4): 939--56.
:::

::: {#ref-Beiko2006 .csl-entry}
Beiko, Robert G, and Nicholas Hamilton. 2006. "[Phylogenetic identification of lateral genetic transfer events.]{.nocase}" *BMC Evolutionary Biology* 6 (January): 15. <https://doi.org/10.1186/1471-2148-6-15>.
:::

::: {#ref-Bininda2004 .csl-entry}
Bininda-Emonds, Olaf R P. 2004. "[The evolution of supertrees.]{.nocase}" *Trends in Ecology & Evolution* 19 (6): 315--22. <https://doi.org/10.1016/j.tree.2004.03.015>.
:::

::: {#ref-Bininda2011 .csl-entry}
---------. 2011. "[Inferring the Tree of Life: chopping a phylogenomic problem down to size?]{.nocase}" *BMC Biology* 9 (Figure 2): 59. <https://doi.org/10.1186/1741-7007-9-59>.
:::

::: {#ref-Bogdanowicz08 .csl-entry}
Bogdanowicz, Damian. 2008. "Comparing Phylogenetic Trees Using a Minimum Weight Perfect Matching." *2008 1st International Conference on Information Technology*, May, 1--4. <https://doi.org/10.1109/INFTECH.2008.4621680>.
:::

::: {#ref-Boussau2013 .csl-entry}
Boussau, Bastien, Gergely J Szöllosi, Laurent Duret, Manolo Gouy, Eric Tannier, and Vincent Daubin. 2013. "Genome-Scale Coestimation of Species and Gene Trees." *Genome Res.* 23 (2): 323--30. <https://doi.org/10.1101/gr.141978.112>.
:::

::: {#ref-Bryant2008 .csl-entry}
Bryant, David, and Mike Steel. 2008. "[Computing the Distribution of a Tree Metric]{.nocase}." *Discrete Mathematics* 6 (3): 420--26. <https://doi.org/10.1109/TCBB.2009.32>.
:::

::: {#ref-Caimo2010 .csl-entry}
Caimo, Alberto, and Nial Friel. 2010. "[Bayesian inference for exponential random graph models]{.nocase}." *Arxiv* 1007.5192 (July): 29. <http://arxiv.org/abs/1007.5192>.
:::

::: {#ref-Caimo2012 .csl-entry}
---------. 2012. "[Bayesian model selection for exponential random graph models]{.nocase}." *Arxiv* 1201.2337 (January): 23. <http://arxiv.org/abs/1201.2337>.
:::

::: {#ref-Chaudhary2010 .csl-entry}
Chaudhary, Ruchi, Mukul S Bansal, André Wehe, David Fernández-Baca, and Oliver Eulenstein. 2010. "[iGTP: a software package for large-scale gene tree parsimony analysis.]{.nocase}" *BMC Bioinformatics* 11 (1): 574. <https://doi.org/10.1186/1471-2105-11-574>.
:::

::: {#ref-Chaudhary2013 .csl-entry}
Chaudhary, Ruchi, J. Gordon Burleigh, and David Fernández-Baca. 2013. "[Inferring species trees from incongruent multi-copy gene trees using the Robinson-Foulds distance]{.nocase}." *Algorithms for Molecular Biology* 8 (1): 28. <https://doi.org/10.1186/1748-7188-8-28>.
:::

::: {#ref-Chung2013 .csl-entry}
Chung, Yujin, Nicole T. Perna, and Cecile Ane. 2013. "Computing the Joint Distribution of Tree Shape and Tree Distance for Gene Tree Inference and Recombination Detection." *IEEE/ACM Transactions on Computational Biology and Bioinformatics*. <https://doi.org/10.1109/TCBB.2013.109>.
:::

::: {#ref-Cotton2009 .csl-entry}
Cotton, James A, and M. Wilkinson. 2009. "[Supertrees join the mainstream of phylogenetics]{.nocase}." *Trends in Ecology & Evolution* 24 (1): 1--3. <http://linkinghub.elsevier.com/retrieve/pii/S0169534708003261>.
:::

::: {#ref-Day1985 .csl-entry}
Day, WHE. 1985. "[Optimal algorithms for comparing trees with labeled leaves]{.nocase}." *Journal of Classification* 2: 7--28. <http://www.springerlink.com/index/q5906x80g44p44k8.pdf>.
:::

::: {#ref-Degiorgio2010 .csl-entry}
Degiorgio, Michael, and James H Degnan. 2010. "[Fast and consistent estimation of species trees using supermatrix rooted triples.]{.nocase}" *Molecular Biology and Evolution* 27 (3): 552--69. <http://www.ncbi.nlm.nih.gov/pubmed/19833741>.
:::

::: {#ref-Doyon2011 .csl-entry}
Doyon, Jean-Philippe, Celine Scornavacca, K Yu Gorbunov, Gergely J Szöllősi, Vincent Ranwez, and Vincent Berry. 2011. "[An efficient algorithm for gene/species trees parsimonious reconciliation with losses, duplications and transfers]{.nocase}." *Comparative Genomics*, no. October: 93--108. <http://www.springerlink.com/index/D8H78633P3N31268.pdf>.
:::

::: {#ref-Drummond2012 .csl-entry}
Drummond, Alexei J., Marc A. Suchard, Dong Xie, and Andrew Rambaut. 2012. "[Bayesian phylogenetics with BEAUti and the BEAST 1.7.]{.nocase}" *Molecular Biology and Evolution* 29 (8): 1969--73. <https://doi.org/10.1093/molbev/mss075>.
:::

::: {#ref-Dunn2008 .csl-entry}
Dunn, Casey W, Andreas Hejnol, D. Q. Matus, K. Pang, W. E. Browne, S. A. Smith, E. Seaver, et al. 2008. "[Broad phylogenomic sampling improves resolution of the animal tree of life]{.nocase}." *Nature* 452 (7188): 745--49. <http://www.nature.com/nature/journal/vaop/ncurrent/abs/nature06614.html>.
:::

::: {#ref-Edwards2007 .csl-entry}
Edwards, Scott V, Liang Liu, and Dennis K Pearl. 2007. "[High-resolution species trees without concatenation.]{.nocase}" *Proceedings of the National Academy of Sciences of the United States of America* 104 (14): 5936--41. <https://doi.org/10.1073/pnas.0607004104>.
:::

::: {#ref-Feng2003 .csl-entry}
Feng, Xizhou, Duncan A. Buell, John R. Rose, and Peter J. Waddell. 2003. "Parallel Algorithms for Bayesian Phylogenetic Inference." *Journal of Parallel and Distributed Computing* 63 (7-8): 707--18. <https://doi.org/10.1016/S0743-7315(03)00079-0>.
:::

::: {#ref-Feng2006 .csl-entry}
Feng, Xizhou, Kirk W. Cameron, and Duncan A. Buell. 2006. "PBPI: A High Performance Implementation of Bayesian Phylogenetic Inference." *Proceedings of the 2006 ACM/IEEE Conference on Supercomputing*. <https://doi.org/10.1145/1188455.1188535>.
:::

::: {#ref-Fletcher2009 .csl-entry}
Fletcher, William, and Ziheng Yang. 2009. "[INDELible: a flexible simulator of biological sequence evolution.]{.nocase}" *Molecular Biology and Evolution* 26 (8): 1879--88. <https://doi.org/10.1093/molbev/msp098>.
:::

::: {#ref-Pawe2012 .csl-entry}
Górecki, Pawe, and Oliver Eulenstein. 2012. "[A Robinson-Foulds Measure to Compare Unrooted Trees with Rooted Trees]{.nocase}." *Bioinformatics Research and Applications* 7292: 115--26. [https://doi.org/10.1007/978-3-642-30191-9\\\_12](https://doi.org/10.1007/978-3-642-30191-9\_12).
:::

::: {#ref-Guigo1996 .csl-entry}
Guigó, R, I Muchnik, and T F Smith. 1996. "[Reconstruction of ancient molecular phylogeny.]{.nocase}" *Molecular Phylogenetics and Evolution* 6 (2): 189--213. <https://doi.org/10.1006/mpev.1996.0071>.
:::

::: {#ref-Heled2010 .csl-entry}
Heled, Joseph, and Alexei J Drummond. 2010. "[Bayesian inference of species trees from multilocus data.]{.nocase}" *Molecular Biology and Evolution* 27 (3): 570--80. <https://doi.org/10.1093/molbev/msp274>.
:::

::: {#ref-Helmkamp2012 .csl-entry}
Helmkamp, Laura J, Ethan M Jewett, and Noah A Rosenberg. 2012. "Improvements to a Class of Distance Matrix Methods for Inferring Species Trees from Gene Trees." *J. Comput. Biol.* 19 (6): 632--49. <https://doi.org/10.1089/cmb.2012.0042>.
:::

::: {#ref-Holton2010 .csl-entry}
Holton, Thérèse a, and Davide Pisani. 2010. "[Deep genomic-scale analyses of the metazoa reject Coelomata: evidence from single- and multigene families analyzed under a supertree and supermatrix paradigm.]{.nocase}" *Genome Biology and Evolution* 2 (January): 310--24. <https://doi.org/10.1093/gbe/evq016>.
:::

::: {#ref-Huggins2011 .csl-entry}
Huggins, P M, W Li, D Haws, T Friedrich, J Liu, and R Yoshida. 2011. "[Bayes estimators for phylogenetic reconstruction.]{.nocase}" *Systematic Biology* 60 (4): 528--40. <https://doi.org/10.1093/sysbio/syr021>.
:::

::: {#ref-Katz2012 .csl-entry}
Katz, L. A., J. R. Grant, L. W. Parfrey, and J. G. Burleigh. 2012. "Turning the Crown Upside down: Gene Tree Parsimony Roots the Eukaryotic Tree of Life." *Systematic Biology* 61 (4): 653--60. <https://doi.org/10.1093/sysbio/sys026>.
:::

::: {#ref-Koonin2008 .csl-entry}
Koonin, Eugene V, and Yuri I Wolf. 2008. "[Genomics of bacteria and archaea: the emerging dynamic view of the prokaryotic world.]{.nocase}" *Nucleic Acids Research* 36 (21): 6688--6719. <https://doi.org/10.1093/nar/gkn668>.
:::

::: {#ref-Kubatko2009b .csl-entry}
Kubatko, Laura Salter, Bryan C Carstens, and L Lacey Knowles. 2009. "[STEM: species tree estimation using maximum likelihood for gene trees under coalescence.]{.nocase}" *Bioinformatics* 25 (7): 971--73. <https://doi.org/10.1093/bioinformatics/btp079>.
:::

::: {#ref-Lang2013 .csl-entry}
Lang, Jenna Morgan, Aaron E. Darling, and Jonathan A. Eisen. 2013. "Phylogeny of Bacterial and Archaeal Genomes Using Conserved Genes: Supertrees and Supermatrices." Edited by Paul J.Editor Planet. *PLoS ONE* 8 (4): e62510. <https://doi.org/10.1371/journal.pone.0062510>.
:::

::: {#ref-Lanier2014 .csl-entry}
Lanier, Hayley C., Huateng Huang, and L. Lacey Knowles. 2014. "How Low Can You Go? The Effects of Mutation Rate on the Accuracy of Species-Tree Estimation." *Molecular Phylogenetics and Evolution* 70 (September): 112--19. <https://doi.org/10.1016/j.ympev.2013.09.006>.
:::

::: {#ref-Larget1999 .csl-entry}
Larget, B, and DL Simon. 1999. "[Markov chain Monte Carlo algorithms for the Bayesian analysis of phylogenetic trees]{.nocase}." *Molecular Biology and Evolution* 16 (6): 750--59. <http://www.renyi.hu/conferences/bp2008/Slides/larget1-budapest.pdf>.
:::

::: {#ref-Lartillot2004 .csl-entry}
Lartillot, Nicolas, and Hervé Philippe. 2004. "[A Bayesian mixture model for across-site heterogeneities in the amino-acid replacement process.]{.nocase}" *Molecular Biology and Evolution* 21 (6): 1095--1109. <https://doi.org/10.1093/molbev/msh112>.
:::

::: {#ref-Leache2011 .csl-entry}
Leaché, Adam D, and Bruce Rannala. 2011. "[The accuracy of species tree estimation under simulation: a comparison of methods.]{.nocase}" *Systematic Biology* 60 (2): 126--37. <https://doi.org/10.1093/sysbio/syq073>.
:::

::: {#ref-Liang2010 .csl-entry}
Liang, Faming. 2010. "[A double Metropolis-Hastings sampler for spatial models with intractable normalizing constants]{.nocase}." *Journal of Statistical Computation and Simulation* 80 (9): 1007--22. <https://doi.org/10.1080/00949650902882162>.
:::

::: {#ref-Liu2000 .csl-entry}
Liu, JS, F Liang, and WH Wong. 2000. "[The multiple-try method and local optimization in Metropolis sampling]{.nocase}." *Journal of the American Statistical Association* 95 (449): 121. <http://www.tandfonline.com/doi/abs/10.1080/01621459.2000.10473908>.
:::

::: {#ref-Liu2001 .csl-entry}
Liu, Jun S. 2001. *Monte Carlo Strategies in Scientific Computing*. Springer.
:::

::: {#ref-Liu2010phybase .csl-entry}
Liu, Liang, and Lili Yu. 2010. "[Phybase: an R package for species tree analysis]{.nocase}." *Bioinformatics* 26 (7): 962--63. <https://doi.org/10.1093/bioinformatics/btq062>.
:::

::: {#ref-Liu2009 .csl-entry}
Liu, Liang, Lili Yu, Laura Salter Kubatko, Dennis K Pearl, and Scott V Edwards. 2009. "[Coalescent methods for estimating phylogenetic trees.]{.nocase}" *Molecular Phylogenetics and Evolution* 53 (1): 320--28. <https://doi.org/10.1016/j.ympev.2009.05.033>.
:::

::: {#ref-Liu2009e .csl-entry}
Liu, L., Lili Yu, Dennis K Pearl, and Scott V Edwards. 2009. "[Estimating Species Phylogenies Using Coalescence Times among Sequences]{.nocase}." *Systematic Biology* 58 (5): 468--77. <https://doi.org/10.1093/sysbio/syp031>.
:::

::: {#ref-Maddison1997 .csl-entry}
Maddison, Wayne P. 1997. "[Gene Trees in Species Trees]{.nocase}." *Systematic Biology* 46 (3): 523--36. <https://doi.org/10.1093/sysbio/46.3.523>.
:::

::: {#ref-Maddison2006 .csl-entry}
Maddison, Wayne P, and L Lacey Knowles. 2006. "[Inferring phylogeny despite incomplete lineage sorting.]{.nocase}" *Systematic Biology* 55 (1): 21--30. <https://doi.org/10.1080/10635150500354928>.
:::

::: {#ref-Mallo2013 .csl-entry}
Mallo, Diego, Leonardo de Oliveira Martins, and David Posada. 2014. "SimPhy: Comprehensive Simulation of Gene, Locus and Species Trees at the Genome-Wide Level." *Unpublished*. <https://bitbucket.org/dmallo/stidsim>.
:::

::: {#ref-Medina2011 .csl-entry}
Medina, Edgar M, Gary W Jones, and David A Fitzpatrick. 2011. "[Reconstructing the fungal tree of life using phylogenomics and a preliminary investigation of the distribution of yeast prion-like proteins in the fungal kingdom.]{.nocase}" *Journal of Molecular Evolution* 73 (3-4): 116--33. <https://doi.org/10.1007/s00239-011-9461-4>.
:::

::: {#ref-Meng2009 .csl-entry}
Meng, Chen, and Laura Salter Kubatko. 2009. "[Detecting hybrid speciation in the presence of incomplete lineage sorting using gene tree incongruence: a model.]{.nocase}" *Theoretical Population Biology* 75 (1): 35--45. <https://doi.org/10.1016/j.tpb.2008.10.004>.
:::

::: {#ref-Mossel2008 .csl-entry}
Mossel, Elchanan, and Sebastien Roch. 2008. "[Incomplete lineage sorting: consistent phylogeny estimation from multiple loci.]{.nocase}" *IEEE/ACM Transactions on Computational Biology and Bioinformatics / IEEE, ACM* 7 (1): 166--71. <https://doi.org/10.1109/TCBB.2008.66>.
:::

::: {#ref-Murray2012 .csl-entry}
Murray, Iain, Zoubin Ghahramani, and David MacKay. 2012. "[MCMC for doubly-intractable distributions]{.nocase}." *Arxiv* 1206.6848 (June). <http://arxiv.org/abs/1206.6848>.
:::

::: {#ref-Nye2006 .csl-entry}
Nye, Tom M W, Pietro Liò, and Walter R Gilks. 2006. "[A novel algorithm and web-based tool for comparing two alternative phylogenetic trees.]{.nocase}" *Bioinformatics* 22 (1): 117--19. <https://doi.org/10.1093/bioinformatics/bti720>.
:::

::: {#ref-Martins10 .csl-entry}
Oliveira Martins, Leonardo de, and Hirohisa Kishino. 2010. "Distribution of Distances Between Topologies and Its Effect on Detection of Phylogenetic Recombination." *Annals of the Institute of Statistical Mathematics* 62 (1): 145--59. <https://doi.org/10.1007/s10463-009-0259-8>.
:::

::: {#ref-Martins08 .csl-entry}
Oliveira Martins, Leonardo de, Elcio Leal, and Hirohisa Kishino. 2008. "Phylogenetic Detection of Recombination with a Bayesian Prior on the Distance Between Trees." *PLoS ONE* 3 (7): e2651. <https://doi.org/10.1371/journal.pone.0002651>.
:::

::: {#ref-Page1994 .csl-entry}
Page, RDM. 1994. "Maps Between Trees and Cladistic Analysis of Historical Associations Among Genes, Organisms and Areas." *Systematic Biology* 43 (1): 58--77. <http://sysbio.oxfordjournals.org/content/43/1/58.short>.
:::

::: {#ref-Pandolfi2010 .csl-entry}
Pandolfi, Silvia, Francesco Bartolucci, and PN Friel. 2010. "[A generalization of the Multiple-try Metropolis algorithm for Bayesian estimation and model selection]{.nocase}." *AISTATS* 9: 581--88. <http://eprints.pascal-network.org/archive/00007723/>.
:::

::: {#ref-Ragan1992 .csl-entry}
Ragan, M. A. 1992. "Phylogenetic Inference Based on Matrix Representation of Trees." *Mol Phylogenet Evol* 1 (1): 53--58. <http://www.ncbi.nlm.nih.gov/pubmed/1342924>.
:::

::: {#ref-Rannala2003 .csl-entry}
Rannala, Bruce, and Ziheng Yang. 2003. "[Bayes Estimation of Species Divergence Times and Ancestral Population Sizes Using DNA Sequences from Multiple Loci.](https://www.ncbi.nlm.nih.gov/pubmed/12930768)" *Genetics* 164 (4): 1645--56.
:::

::: {#ref-Rasmussen2011 .csl-entry}
Rasmussen, Matthew D, and Manolis Kellis. 2011. "[A Bayesian approach for fast and accurate gene tree reconstruction.]{.nocase}" *Molecular Biology and Evolution* 28 (1): 273--90. <https://doi.org/10.1093/molbev/msq189>.
:::

::: {#ref-Rasmussen2012 .csl-entry}
---------. 2012. "[Unified modeling of gene duplication, loss, and coalescence using a locus tree]{.nocase}." *Genome Research* 22 (4): 755--65. <https://doi.org/10.1101/gr.123901.111>.
:::

::: {#ref-Robinson1981 .csl-entry}
Robinson, DF, and LR Foulds. 1981. "[Comparison of phylogenetic trees]{.nocase}." *Mathematical Biosciences* 53: 131--47. <http://www.sciencedirect.com/science/article/pii/0025556481900432>.
:::

::: {#ref-Romiguier2013 .csl-entry}
Romiguier, Jonathan, Vincent Ranwez, Frédéric Delsuc, Nicolas Galtier, and Emmanuel J P Douzery. 2013. "Less Is More in Mammalian Phylogenomics: AT-Rich Genes Minimize Tree Conflicts and Unravel the Root of Placental Mammals." *Mol. Biol. Evol.* 30 (9): 2134--44. <https://doi.org/10.1093/molbev/mst116>.
:::

::: {#ref-Ronquist2004 .csl-entry}
Ronquist, Fredrik, John P. Huelsenbeck, and Tom Britton. 2004. "Bayesian Supertrees." In *Phylogenetic Supertrees*, edited by Olaf R. P. Bininda-Emonds, 193--224. Springer.
:::

::: {#ref-Ronquist2012 .csl-entry}
Ronquist, Fredrik, Maxim Teslenko, Paul van der Mark, Daniel L Ayres, Aaron Darling, Sebastian Höhna, Bret Larget, Liang Liu, Marc A Suchard, and John P Huelsenbeck. 2012. "[MrBayes 3.2: Efficient Bayesian Phylogenetic Inference and Model Choice Across a Large Model Space.]{.nocase}" *Systematic Biology*, February, 1--14. <https://doi.org/10.1093/sysbio/sys029>.
:::

::: {#ref-Rubenthaler09 .csl-entry}
Rubenthaler, Sylvain, Tobias Rydén, and Magnus Wiktorsson. 2009. "Fast Simulated Annealing in Rd with an Application to Maximum Likelihood Estimation in State-Space Models." *Stochastic Processes and Their Applications* 119 (6): 1912--31. <https://doi.org/10.1016/j.spa.2008.09.007>.
:::

::: {#ref-Salichos13 .csl-entry}
Salichos, Leonidas, and Antonis Rokas. 2013. "Inferring Ancient Divergences Requires Genes with Strong Phylogenetic Signals." *Nature* 497 (7449): 327--31. <https://doi.org/10.1038/nature12130>.
:::

::: {#ref-Sanderson2011 .csl-entry}
Sanderson, Michael J, Michelle M McMahon, and Mike Steel. 2011. "[Terraces in phylogenetic tree space.]{.nocase}" *Science (New York, N.Y.)* 333 (6041): 448--50. <https://doi.org/10.1126/science.1206357>.
:::

::: {#ref-Schreiber2013 .csl-entry}
Schreiber, F., M. Patricio, M. Muffato, M. Pignatelli, and A. Bateman. 2014. "TreeFam V9: A New Website, More Species and Orthology-on-the-Fly." *Nucleic Acids Research* 42 (D1): D922--25. <https://doi.org/10.1093/nar/gkt1055>.
:::

::: {#ref-Sjostrand2012 .csl-entry}
Sjostrand, J., B. Sennblad, L. Arvestad, and J. Lagergren. 2012. "[DLRS: gene tree evolution in light of a species tree]{.nocase}." *Bioinformatics* 28 (22): 2994--95. <https://doi.org/10.1093/bioinformatics/bts548>.
:::

::: {#ref-smith1992 .csl-entry}
Smith, Adrian FM, and Alan E Gelfand. 1992. "Bayesian Statistics Without Tears: A Sampling--Resampling Perspective." *The American Statistician* 46 (2): 84--88. <http://www.jstor.org/stable/2684170>.
:::

::: {#ref-Stark2007 .csl-entry}
Stark, Alexander, Michael F. Lin, Pouya Kheradpour, Jakob S. Pedersen, Leopold Parts, Joseph W. Carlson, Madeline A. Crosby, et al. 2007. "Discovery of Functional Elements in 12 Drosophila Genomes Using Evolutionary Signatures." *Nature* 450 (7167): 219--32. <https://doi.org/10.1038/nature06340>.
:::

::: {#ref-Steel1993 .csl-entry}
Steel, M. A., and David Penny. 1993. "[Distributions of Tree Comparison Metrics--Some New Results]{.nocase}." *Systematic Biology* 42 (2): 126--41. <https://doi.org/10.1093/sysbio/42.2.126>.
:::

::: {#ref-Steel2008 .csl-entry}
Steel, Mike, and Allen Rodrigo. 2008. "[Maximum likelihood supertrees.]{.nocase}" *Systematic Biology* 57 (2): 243--50. <https://doi.org/10.1080/10635150802033014>.
:::

::: {#ref-Sukumaran2010 .csl-entry}
Sukumaran, Jeet, and Mark T Holder. 2010. "DendroPy: A Python Library for Phylogenetic Computing." *Bioinformatics* 26 (12): 1569--71. <https://doi.org/10.1093/bioinformatics/btq228>.
:::

::: {#ref-Szollosi2013 .csl-entry}
Szollosi, G. J., W. Rosikiewicz, B. Boussau, E. Tannier, and V. Daubin. 2013. "Efficient Exploration of the Space of Reconciled Gene Trees." *Systematic Biology* 62 (6): 901--12. <https://doi.org/10.1093/sysbio/syt054>.
:::

::: {#ref-Than2009 .csl-entry}
Than, Cuong, and Luay Nakhleh. 2009. "[Species tree inference by minimizing deep coalescences.]{.nocase}" *PLoS Computational Biology* 5 (9): e1000501. <https://doi.org/10.1371/journal.pcbi.1000501>.
:::

::: {#ref-Tofigh2011 .csl-entry}
Tofigh, Ali, Michael Hallett, and Jens Lagergren. 2011. "[Simultaneous identification of duplications and lateral gene transfers.]{.nocase}" *IEEE/ACM Transactions on Computational Biology and Bioinformatics (TCBB)* 8 (2): 517--35. <https://doi.org/10.1109/TCBB.2010.14>.
:::

::: {#ref-Waterman1978 .csl-entry}
Waterman, M. S., and T. F. Smith. 1978. "On the Similarity of Dendrograms." *Journal of Theoretical Biology* 73 (4): 789--800. <https://doi.org/10.1016/0022-5193(78)90137-6>.
:::

::: {#ref-Whidden2013 .csl-entry}
Whidden, Chris, Norbert Zeh, and Robert G. Beiko. 2013. "Supertrees Based on the Subtree Prune-and-Regraft Distance." *PeerJ PrePrints* 1 (May): 18. <https://doi.org/10.7287/peerj.preprints.18v1>.
:::

::: {#ref-Williams2012 .csl-entry}
Williams, Tom A, Peter G Foster, Tom M W Nye, Cymon J Cox, and T Martin Embley. 2012. "[A congruent phylogenomic signal places eukaryotes within the Archaea.]{.nocase}" *Proceedings of the The Royal Society / Biological Sciences* 279 (1749): 4870--79. <https://doi.org/10.1098/rspb.2012.1795>.
:::

::: {#ref-Wu2012b .csl-entry}
Wu, Martin, Sourav Chatterji, and Jonathan a. Eisen. 2012. "[Accounting For Alignment Uncertainty in Phylogenomics]{.nocase}." Edited by Marco Salemi. *PLoS ONE* 7 (1): e30288. <https://doi.org/10.1371/journal.pone.0030288>.
:::

::: {#ref-Wu2011 .csl-entry}
Wu, Taoyang, and Louxin Zhang. 2011. "[Structural properties of the reconciliation space and their applications in enumerating nearly-optimal reconciliations between a gene tree and a species tree.]{.nocase}" *BMC Bioinformatics* 12 Suppl 9 (Suppl 9): S7. <https://doi.org/10.1186/1471-2105-12-S9-S7>.
:::

::: {#ref-Wu2013 .csl-entry}
Wu, Y.-C., M. D. Rasmussen, M. S. Bansal, and M. Kellis. 2013. "Most Parsimonious Reconciliation in the Presence of Gene Duplication, Loss, and Deep Coalescence Using Labeled Coalescent Trees." *Genome Research*, December. <https://doi.org/10.1101/gr.161968.113>.
:::

::: {#ref-Yang2007 .csl-entry}
Yang, Ziheng. 2007. "[PAML 4: phylogenetic analysis by maximum likelihood.]{.nocase}" *Molecular Biology and Evolution* 24 (8): 1586--91. <https://doi.org/10.1093/molbev/msm088>.
:::

::: {#ref-Yang1997 .csl-entry}
Yang, Ziheng, and Bruce Rannala. 1997. "[Bayesian phylogenetic inference using DNA sequences: a Markov Chain Monte Carlo Method.]{.nocase}" *Molecular Biology and Evolution* 14 (7): 717--24. <http://www.ncbi.nlm.nih.gov/pubmed/9214744>.
:::

::: {#ref-Yu2011 .csl-entry}
Yu, Yun, Cuong Than, James H Degnan, and Luay Nakhleh. 2011. "[Coalescent histories on phylogenetic networks and detection of hybridization despite incomplete lineage sorting.]{.nocase}" *Systematic Biology* 60 (2): 138--49. <https://doi.org/10.1093/sysbio/syq084>.
:::

::: {#ref-Zhang1997 .csl-entry}
Zhang, L. 1997. "[On a Mirkin-Muchnik-Smith conjecture for comparing molecular phylogenies]{.nocase}." *Journal of Computational Biology* 4 (2): 177--87. <http://online.liebertpub.com/doi/abs/10.1089/cmb.1997.4.177>.
:::

::: {#ref-Zmasek2001 .csl-entry}
Zmasek, CM, and SR Eddy. 2001. "[A simple algorithm to infer gene duplication and speciation events on a gene tree]{.nocase}." *Bioinformatics* 17 (9): 821--28. <http://bioinformatics.oxfordjournals.org/content/17/9/821.short>.
:::
:::

# Exact source archive
The fenced blocks below are exact UTF-8 source text, not additional article sections. SHA-256 hashes refer to the original file bytes.

## Source: paper001_guenomu/20140204.manuscript.tex
SHA-256: `d51f1c0d6d843a09f6861c4d226cc52d23ab87dfd0264a0aacf8142ac2505491`

````latex
\documentclass[pdftex,12pt,a4paper,openany,fleqn]{article}
\usepackage{indentfirst}
\usepackage[dvips]{graphicx}
\usepackage[dvips,usenames,dvipsnames]{color}
%\usepackage[sort&compress,numbers]{natbib} % bibliographic style (sorted numbers)
\usepackage{natbib} % bibliographic style (sorted numbers)
%\usepackage[sortlocale=en_US,backend=biber,style=draft]{biblatex}
\usepackage{nameref}  % standard references (to chapter, page)
\usepackage{times}    % recommended fonts for pdflatex
\usepackage{mathpazo} % recommended fonts for pdflatex
\usepackage{courier}  % recommended fonts for pdflatex
\usepackage{helvet}   % recommended fonts for pdflatex
\usepackage{cmap}     % to produce searchable PDF
\usepackage{amsmath}
\usepackage{cancel}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{amsthm}
\usepackage{multirow}  % fields that span multiple rows in tabular environment
\usepackage{bm}       % bold math symbols (use as "{\bm \alpha}" or "\bm{\alpha}")
\usepackage{setspace} % use as \begin{spacing}{0.8} to reduce linespacing to 80% 
\usepackage{fancyhdr} % headers and footnotes 
\usepackage{fancyvrb} % \begin{Verbatim} environment
\usepackage{sidecap}  % puts figure/table and its caption side by side
\usepackage[usenames]{color} % colored fonts ("usenames" defines 68 colors  http://en.wikibooks.org/wiki/LaTeX/Colors)

%\usepackage[disable]{todonotes} % notes not showed
\usepackage[draft]{todonotes}  % notes showed (use as \todo{comment} or \todo[inline]{comment} )


\usepackage[sf,small,compact]{titlesec}  % space saving title sectioning
%\usepackage[font=small]{caption, subfig} % caption and subfigure (composite figs) control 
\usepackage[top=0.4in, text={7.2in,10.5in}]{geometry}	% page layout

\usepackage[pdftex, colorlinks, pdfauthor={Leonardo de Oliveira Martins, Diego Mallo, David Posada},
						pdftitle={Bayesian supertree models for genome-wide species tree reconstruction},
						pdfkeywords={},
						pdfsubject={},
						pdfproducer={Latex with hyperref on Debian GNU/Linux},
						pdfcreator={pdflatex}, citecolor=blue, urlcolor=blue, linkcolor=blue]{hyperref}

\pdfinfo{
	/Title      (A Bayesian Supertree Model for Genome-Wide Species Tree Reconstruction)
	/Keywords   ()
	/Author     (Leonardo de Oliveira Martins, Diego Mallo, David Posada)
}

\renewcommand\floatpagefraction{.9} %page fraction with figures/tables 
\renewcommand\topfraction{.9}       %page fraction with figures/tables 
\renewcommand\bottomfraction{.9}    %page fraction with figures/tables 
\renewcommand\textfraction{.1}      %page fraction with text 
\setcounter{totalnumber}{50}  % max figures per page
\setcounter{topnumber}{50}    % max figures per page
\setcounter{bottomnumber}{50} % max figures per page


%title unused, but writelatex needs it to identify doc name
\title{Bayesian supertree models for species tree reconstruction}

% color definitions
\definecolor{blucol}{rgb}{0,0.4,0.7}
\newcommand{\blue}{\textcolor{blucol}}
\newcommand{\green}{\textcolor{green}}
\newcommand{\red}{\textcolor{red}}
\newcommand{\magenta}{\textcolor{magenta}}
\newcommand{\cyan}{\textcolor{cyan}}
\definecolor{orangecol}{rgb}{1,0.5,0}
\newcommand{\orange}{\textcolor{orangecol}}
\newcommand{\deriv}{\;\mathrm{d}}


%\linespread{1.4}
\linespread{1.66} %( for systbiol)  All text should be double-spaced  % with occasional exceptions for tables. 
\raggedright
\setlength{\parindent}{0.5in}
%\bibliographystyle{plainnat}
%\addbibresource{bibtex_probpars.bib}
\bibliographystyle{sysbio}
\bibpunct{(}{)}{;}{a}{}{,}  % this is a citation format command for natbib


\setcounter{secnumdepth}{0}
% SystBiol sections are not numbered and our papers do not have Tables of Contents. We don't  present a list of figures or list of tables, either.
% Any common font is fine. (A common sans-serif font should be used on figures, but figures should be separate from the LaTeX document.)
\pagestyle{empty}
\renewcommand{\section}[1]{%
\bigskip\begin{center}\begin{Large} \normalfont\scshape #1 \medskip \end{Large}\end{center}}
\renewcommand{\subsection}[1]{%
\bigskip\begin{center}\begin{large} \normalfont\itshape #1 \end{large}\end{center}}
\renewcommand{\subsubsection}[1]{%
\vspace{2ex} \noindent \textit{#1.}---}
\renewcommand{\tableofcontents}{}

\begin{document}
\begin{flushright} Version dated: \today  \end{flushright}
\bigskip \noindent RH: BAYESIAN SUPERTREE RECONSTRUCTION 
\bigskip \medskip
\begin{center}

\noindent{\Large \bf A Bayesian Supertree Model for Genome-Wide Species Tree Reconstruction}
\bigskip
% We don't use a special title page; the author information is entered  like any other text.
% FOOTNOTES: We don't allow them in the manuscript, except in tables. Don't include any footnotes in the text.

\noindent {\normalsize \sc Leonardo de Oliveira Martins$^1$, Diego Mallo$^2$, and David Posada$^3$}\\
\noindent {\small \it  Department of Biochemistry, Genetics and Immunology, University of Vigo, Vigo, 36310,  Spain}\\
%Tel +34 986 818 681
%Fax: +34 986 812 556
\end{center}
\medskip
\noindent{\bf Corresponding author:} Leonardo de Oliveira Martins, 
Department of Biochemistry, Genetics and Immunology, University of Vigo, Vigo, 36310,  Spain\\

% Of course the specific format of addresses may vary according to country or other factors. Also, that was just an example email format.
%It's acceptable to add email addresses for authors in addition to the corresponding author. These would be placed after "Country."
\noindent{\bf E-mail:} $^1$leomrtns@uvigo.es, \quad $^2$ dmallo@uvigo.es, \quad $^3$ dposada@uvigo.es\\
\vspace{0.6in}

\subsubsection{Abstract}

Current phylogenomic data sets highlight the need for species tree methods able to deal with several sources of gene
tree/species tree incongruence. At the same time, we need to make most use of all available data. Most species tree
methods deal with single processes of phylogenetic discordance, namely, gene duplication and loss, incomplete lineage
sorting or horizontal gene transfer. In this manuscript we address for the first time the problem of species tree
inference from multilocus, genome-wide data sets in the  presence of gene duplication and loss and incomplete lineage
sorting, therefore without the need to identify orthologs or to use a single individual per species. We do this by
extending the idea of Maximum Likelihood supertrees to a hierarchical Bayesian model where several sources of gene
tree/species tree disagreement can be accounted for in a modular manner.  We implemented this model in a computer
program called {\em guenomu} whose input are with posterior distributions of unrooted gene tree topologies for multiple
gene families, and whose output is the posterior distribution of rooted species tree topologies. We  conducted extensive
simulations under complex phylogenomic models in order to evaluate the performance ouf our approach in comparison with
other species tree approaches able to deal with multilabeled trees. 
Our method ranked best under simulated data sets, in spite of ignoring branch lengths, and performed well on empirical
data sets as well.
Our results show in addition  that under complex simulation
scenarios , gene tree parsimony is also a competitive approach once we consider its speed in contrast to more
sophisticated models.

\noindent (Keywords: supertree, reconciliation, hierarchical Bayesian model, tree distance, phylogenomics)\\
\vspace{1in}

The evolutionary history of a gene family is not an exact representation of the evolution of the species embedding it,
since we must take into account the effects of gene duplication, population factors like incomplete lineage sorting as
well as the possibility of horizontal gene transfer \citep{Maddison1997}. 
In recent years the rapid accumulation of phylogenomic data has highlighted the need for specific methods to infer
species phylogenies.  
Methods of species tree inference based on collection of gene alignments/trees can be broadly classified into
supermatrix, supertree and model-based approaches.
The supermatrix strategy consists in concatenating different genes into a single large alignment, assuming that minor
conflicting phylogenetic signals will cancel out and a common phylogenetic trend -- the species history -- will emerge.
Departures from this assumption must be handled by subjectively removing columns or whole regions \citep{Wu2012b,
Lanier2014}.  This method relies on the correct inference of the orthologous sets, such that the sequences from a single
species can be concatenated into a single 'supergene' \citep{Degiorgio2010}. 

Alternatively, phylogenetic inference can be conducted independently for each gene alignment, and then in a second step
the species tree is inferred from the resulting gene trees \citep{Cotton2009}. In this case the alignments are used
indirectly, through summary statistics like maximum likelihood trees, to build a 'supertree'. Under the supertree
umbrella, some methods try to find the species tree that minimizes its overall disagreement with the gene tree
collection without explicitly  taking into account the biological phenomena behind the incongruence, like the RF
supertree \citep{Bansal2010, Chaudhary2013} and the MRP approaches \citep{Ragan1992,Bininda2004}.  Other supertree
approaches aim for the most parsimonious scenario by assuming that the disagreement between the species tree and each of
the gene trees is due to HGT \citep{Whidden2013}, duplications and losses or deep coalescences \citep{Bansal13}. The
problem of finding the species tree that minimizes the reconciliation cost is called Gene Tree Parsimony (GTP)
\citep{Bansal13}.

Most supertree approaches  neglect the gene tree branch lengths, but there is another class of species tree methods
closely related to the supertree methods which are based on the multispecies coalescent.
These methods try to reconstruct the species tree from a matrix of distances between species, which in turn is built
based on distance matrices from the individual gene trees (reviewed in \cite{Liu2009, Helmkamp2012}).
These methods are based on the observation that the coalescence time of two genes from distinct species always precedes
the speciation event. 
In general they assume that the gene trees are known without error \citep{Kubatko2009b,Liu2010phybase}. Some of
these methods also can take into account hibridization \citep{Meng2009} or HGT \citep{Yu2011}. 

In principle, the most powerful species tree methods are based on an explicit probabilistic modelling of the
incongruence between species and gene trees based on gene alignments. Together with the phylogenetic likelihood that
takes into account the stochastic uncertainty in the gene tree inference, these models also describe the probability of
such gene tree being generated by a given species tree. If we assume that the sequences represent individuals from a
given orthologous region, then this probability is given by the multispecies coalescent, that describes how lineages
coalesce conditioned on a species tree \citep{Rannala2003}. Under this model  Bayesian methods have already been
implemented to estimate the posterior distribution of species trees, together with ancestral population sizes
 and divergence times \citep{Edwards2007, Heled2010, Ronquist2012}.
If, however, we assume that each alignment represents
a set of paralogous gene copies  from a gene family, then we can model the duplications and losses inside each gene
family through a birth-death process \citep{Arvestad2003}. Although theoretically it is possible to estimate the
distribution of species trees while integrating over the particular birth-death histories, current implementations
assume that the species tree is known \citep{Akerborg2009,Rasmussen2011, Sjostrand2012} or try to find it by  maximum
likelihood \citep{Boussau2013}. 
Under the fixed tree assumption, it is furthermore possible to add HGT to the birth-death model \citep{Szollosi2013}.
Fully probabilistic models, although more realistic, tend to be slow. Furthermore, they are usually
limited to particular sources of gene tree / species tree disagreement. Under the multispecies coalescent, it is assumed
that all genes from a same species are orthologous, neglecting the possibility of duplication, while duplication and
loss models assume that sequences mapped to one species are necessarily the product of a duplication
\citep{Rasmussen2012}.
 
Recently a maximum likelihood (ML) supertree approach has been proposed, such that the probabilities are associated to
errors in the gene trees \citep{Steel2008, Cotton2009}. That is, that given the supertree the probability will reflect a
penalty against incongruent input trees, regardless of the evolutionary process responsible for the incongruence. The
model assumes that a true supertree generates trees with error, such that their probabilities decrease exponentially
according to the distance from the original supertree. The incorporation of supertree approaches into a likelihood
framework allows us to compare alternative solutions statistically and to be more explicit about our assumptions,
incorporating model selection and hypothesis testing. The probability $P(T_i\mid\tau,\lambda)$ of supertree $\tau$
generating tree $T_i$ is given by
\begin{equation} 
%        P(T_i\mid\tau,\lambda) = \frac{e^{-d(T_i,\tau)/\lambda}}{Z(\tau,\lambda)}
        P(T_i\mid\tau,\lambda) \propto e^{-d(T_i,\tau)/\lambda}
\end{equation}

where the scale parameter $\lambda$ controls how strongly discordant trees are penalized. The distance $d()$ is usually
a discrete tree-to-tree metric that does not take branch lengths into account, and therefore the original ML supertree
approach is based on tree topologies only. Furthermore, it was always assumed that each input tree contained a subset of
the species, as is commonplace for supertree methods. The particular distance employed  will depend on the problem at
hand, as e.g. the Robinson-Foulds distance will lead to a supertree in agreement  with the  largest number of clusters
(or clades) from the input trees \citep{Bansal2010}. 
Likewise, the subtree prune-regraft (SPR) distance \citep{Beiko2006}
has been used in an exponential model of phylogenetic recombination, where the penalty is against trees in adjacent
alignment segments that can only be reconciled through several recombination events, as estimated by the SPR distance
between them \citep{Martins08, Martins10}.

Probabilistic supertree approaches like ML supetrees have a lot of potential, as they can leverage statistical
sophistication with computational tractability, being able to deal with several processes of incongruence at once
\citep{Cotton2009}. In this work we further extended the ML supertree approach in several notable ways, creating a
hierarchical Bayesian supertree model \citep{Ronquist2004}. Firstly, we use a multivariate description of the
disagreement between each input tree and the supertree, such that several distances can be used at once. Importantly,
among the distances we use jointly are the most parsimonious reconciliation costs under the DL and the ILS models
\citep{Guigo1996,Than2009}.  In our model the input trees represent gene families and the supertree is the species
phylogeny. This allows us to work with more general data sets than supertree methods usually handle, since we can work
with input trees with several leaves representing the same species (even larger, in fact, than the supertree).
Furthermore we incorporate this multivariate distribution of tree-supertree distances into a hierarchical Bayesian
model, such that the uncertainty of the gene trees as well as in the scale parameters $\lambda$ are taken into account.

Our full model can include the phylogenetic likelihood describing how gene trees are supported by the alignments, but
in practice we adopted an importance sampling step such that any Bayesian phylogenetic program can provide the
individual distribution of gene trees.
(That is why we call our method a ’supertree’, in connection to the ML supertree \citep{Steel2008}, despite the full 
model can in principle work with sequence alignments as well). 
The resulting hierarchical Bayesian model accounts then
for the uncertainty within each gene family to estimate the posterior distribution of all compatible species trees.
Unlike the fully probabilistic models described above, our approach is based on the simplifying assumption that
summaries based on the most parsimonious reconciliation scenarios are sufficient to explain the disagreement between
gene and species trees, regardless of its causes. 
Yet we do not rely on known gene trees, as with most supertree approaches. 
Instead we work with marginal gene tree
distributions from  independent Bayesian phylogenetic analyses of multilocus alignments. 
We implemented  this model in a program
called {\em guenomu}, that receives as input a set of unrooted gene tree distributions, one per gene family, as well as
a list of species names, and returns the posterior distribution of rooted species trees as well as a posterior
distribution of gene trees  for each gene family. 


In this manuscript we evaluate our method under simulated scenarios in which gene trees evolve inside species tree with
gene duplications and losses as well as incomplete lineage sorting, and assuming gene tree inference error. Usually,
these different sources of gene tree disagreement are considered separatedly, although unrecognized processes can have a
deleterious effect on the species tree estimation \citep{Rasmussen2012}, since even a single ancestral duplication can
mimic the effect of several deep coalescences.
We circumvent this by favoring species trees that are compatible with
different phenomena at once, as dictated by our multivariate distance distribution. 
Each gene family, on its turn, is allowed to have different contributions in each regard, since we do not expect
the contribution of duplications and losses, for instance, to be uniform along the genome. Our simulations suggest that our
Bayesian hierarchical model is able to reduce gene tree uncertainty and results in reasonable species trees, considering
the complexity of the simulation scenario. Unexpectedly, we also found that GTP methods can be quite accurate as well,
while coalescent summary statistics methods performed badly under the scenarios considered.

\section{Materials and Methods}
\subsection{Definitions}
We define {\it gene family} $i$ as the set of homologous sequences that compose an alignment $D_i$, which can be
comprised of paralogs and orthologs belonging to the same or distinct individuals from a diverged population, that we
will call a {\it species}. The only information needed is a mapping between each member of this gene family and the
species which it belongs to (species are thus defined as the entities whose phylogenies we are interested in).
Importantly, we assume that we do not know beforehand which members of a gene family are related through a duplication
and which diverged trhough coalescences. That is, we  map the leaves from the gene tree  directly to the species tree.
We will then call {\it gene tree} any phylogenetic tree connecting all members of a gene family, whose corresponding
random variable for gene family $i$ will be represented by $G_i$. We are particularly interested in the posterior
distribution of gene trees as e.g. estimated by Bayesian phylogenetic inference programs like MrBayes or PhyloBayes
\citep{Ronquist2012, Lartillot2004}. These will be our {\it input gene tree distributions}, one per gene family, and for
each input distribution of gene trees we will have an associated {\it posterior distribution}, in a manner analogous to
the BUCKy model \citep{Ane2007}. Although we will mostly use the term {\it tree}, the input and output distributions
contain only the tree topologies. That is, in the current implementation of our model we do not use the branch length
information, unless otherwise stated. We will resort to the terms {\it phylogeny} and {\it topology} when we want to be
explicit about the presence or absence, respectively, of branch lengths.


\subsection{Bayesian hierarchical model for species trees}
We have devised a Bayesian hierarchical model in which the posterior distribution of the species tree $S$ and all other
parameters $\Theta$ given a set of $N$ gene family alignments $\bm{D}=(D_1,...,D_N)$ can be described by 
\begin{equation}\label{eq:postprob}
P(S,{\bm \Theta}\mid\bm{D})\propto
P(\bm{D},\bm{\theta}\mid\bm{G}) P(\bm{G}\mid\bm{\lambda},S) P(\bm{\lambda}\mid\bm{\lambda_0}) P(\bm{\lambda_0}) P(S)
\end{equation}
where $G=(G_1,\hdots, G_N)$ are the gene trees and $\bm{\lambda}=(\lambda_{ij})$ is the vector with distance scale parameters
$j$ for gene family $i$, which depend on the hyperhyperparameters $\bm{\lambda_0}$. The term $P(\bm{D},\bm{\theta}\mid \bm{G})$  
is the joint probability of alignments $D$ and parameters $\bm{\theta}$ related to the substitution model -- as we will
see this calculation is delegated to the user, so our model works equally for DNA, codon and protein
alignments.  
All species trees $S$ are equally probable a priori, leading to the uniform prior $P(S)=|\mathcal{S}|^{-1}$ over
all possible rooted trees $S\in\mathcal{S}$.
Equation \ref{eq:postprob} can be further
decomposed into 
\begin{eqnarray*}
  P(S,{\bm \Theta}\mid\bm{D})&\propto&  
        P(\bm{\lambda_0})P(S) \times \prod_{i=1}^N P(D_i,{\bm \theta_i}\mid G_i)
        P(G_i\mid\lambda_{i\cdot},S)P(\lambda_{i\cdot}\mid\bm{\lambda_0})
\end{eqnarray*}
where we see that the hyperhyperparameters $\bm{\lambda_0}$ and the species tree $S$ are shared among all gene families, while other
parameters are allowed to vary between families. 

\cite{Steel2008} assumed that the probability $P(G_i\mid S)$ of a gene tree $G_i$  given the species tree $S$
followed an exponential distribution that depended solely on the dissimilarity between $G_i$ and $S$, which can be represented by a
distance $d(G_I,S)$ and a scale parameter $\lambda_i$. We generalized this distribution and assume that this probability
depends on the disagreement between $G_i$ and $S$ with respect to several distances
$d(G_i,S)=(d_1(G_i,S),\hdots,d_J(G_i,S))$ where each measure $d_j$ $(j=1,...,J)$ represents a distinct distance. 
Therefore the distribution can be written as
\begin{equation}
        P(G_i\mid S) = 
        \frac{e^{-\sum_{j=1}^J d_j(G_i,S)/(m_{ij}\lambda_{ij})}}{\sum_{g\in\mathcal{G}_i} e^{-\sum_{j=1}^J d_j(g,S)/(m_{ij}\lambda_{ij})}} =
        \frac{e^{-\sum_{j=1}^J d_j(G_i,S)/(m_{ij}\lambda_{ij})}}{Z_i(S,\lambda_{i\cdot})} 
\end{equation}

where $Z_i(S,\lambda_{i\cdot})$ is the {\it normalization constant} (also called partition function; see below) 
over the set $G_i$  of all possible topologies for gene family $i$, and the $m_{ij}$ are fixed
constants.

Each gene family $i$ has a vector $\lambda_{i\cdot}=(\lambda_{i1},\hdots,\lambda_{iJ})$ where each parameter 
$\lambda_{ij}$ is associated to a distance $d_j$. We have explored two groups of distances, which we will call {\it
reconciliation distances} and {\it non-parametric distances} (although they are not in general proper metrics since they
do not satisfy the symmetry condition). 
Our model allows for different parameterizations depending on which distances we want to consider for a specific
analysis.
The reconciliation distances are based on the most parsimonious reconciliations between the rooted species tree and gene
trees \citep{Page1994,Maddison1997}, which can be calculated using the least common ancestor (LCA) mapping between
each node $g$ of the gene tree and a corresponding node $M(g)$ of the species tree \citep{Zmasek2001, Bansal13}. This mapping
can thus be used in two independent ways: to find the minimum number of duplications and losses needed to explain the
disagreement \citep{Guigo1996,Zhang1997}, or to find the minimum number of deep coalescences necessary to make the gene
tree be congruent with the species tree \citep{Than2009,Wu2011}.

Under the duplication-loss reconciliation model and given the mapping $M(\cdot)$ between nodes of the gene tree $G_i$
rooted at $r$, denoted $G_i^r$, and the rooted species tree $S$, a gene tree node $g$ will represent a duplication if
$M(g)=M(g_1)$ or $M(g)=M(g_2)$, where $g_1$ and $g_2$ are the two children of node $g$. 
The number of duplications $\mathrm{DUPS}(G_i^r,S)$ will be the number of such nodes, and the number of losses
$\mathrm{LOSS}(G_i^r,S)$ is calculated in postorder based on the LCA mapping and the duplication nodes. 
Likewise, under the reconciliation coalescent model we assume all gene tree nodes are coalescences,
and we associate to each node $s$ on the species tree the number of extra lineages, which is the number of branches from
the gene tree passing through $s$ minus one. 
The number of deep coalescences $\mathrm{ILS}(G_i^r,S)$ -- also called incomplete lineage sorting events -- is the sum
of extra lineages over all species tree nodes. Since we are using the LCA mapping these numbers will be minimal over all
possible reconciliations -- for more details, see \cite{Wu2011, Bansal13}. 
Since we do not know the root location of the gene trees, we apply the LCA mapping and calculate the minimum costs over all
possible roots $r$ to define the optimal root location. From the duplication-loss model we will then have $d_1(G_i,S)=\mathrm{DUPS}(G_i^{r^*},S)$ and
$d_2(G_i,S)=\mathrm{LOSS}(G_i^{r^*},S)$ which will be the smallest amongst all $r$, and for for the coalescent model we
will have $d_3(G_i,S)=\mathrm{ILS}(G_i^{r^{**}},S)$, where the optimal root $r^{**}$ might be the same as $r^{*}$ or not.

As for the non-parametric distances, they can be any estimate of disagreement between rooted or unrooted trees, of which
the most popular example is the Robinson-Foulds (RF) (or symmetric) distance \citep{Robinson1981}. 
Other examples include the SPR distance \citep{Beiko2006} and minimum weight matchings between tree branches \citep{Nye2006, Bogdanowicz08}.
They are non-parametric in the
sense that they do not try to model the biological reason for the disagreement, only the outcome. The RF distance can be
calculated only when each leaf from one tree is mapped to at most one tree from the other, restricting therefore the
gene families that can be considered. 
However very recently a generalization of the RF distance called $mulRF$ was introduced that relaxes this
constraint, allowing for one of the trees to have several leaves with the same label \citep{Chaudhary2013}. 
Such multi-labelled tree, or mul-tree, is here the gene tree. 
Although the original RF distance is between unrooted topologies, the RF distance can be used to compare a rooted tree
to another one, rooted or unrooted \citep{Day1985,Pawe2012}.
Both the $RF$ and the $mulRF$ distances are calculated as the sum of bipartitions present in one tree but not in the
other, but for the $mulRF$ distance the species tree is extended such that all leaves representing a species labelled
more than once in the gene tree are replaced by a multifurcation with same number of leaves as in the gene tree
\citep{Chaudhary2013}. 
Therefore we have that $d_4(G_i,S)=\mathrm{mulRF}(G_i,S)$, noting that this value is the same for any rooting of $G_i$
and S since we use the unrooted version of the trees.
Like the reconciliation distances, the non-parametric ones ignore branch lengths, since they measure the disagreement
between the topologies only.

Our hierarchical model is completed by specifying the priors for the distance scale parameters.
For each distance $j$, the scale parameter $\lambda_{ij}$ of gene family $i$ will follow an exponential hyperprior distribution
\begin{equation}
        P(\lambda_{ij}\mid \bm\lambda_0)=\frac{e^{-\lambda_{ij}/\lambda_{0j}}}{\lambda_{0j}}
\end{equation}
where the scale parameter $\lambda_{0j}$ is the $j$-th element of the (hyperhyperparameter) vector $\bm\lambda_0$,
shared among gene families.
As mentioned in \cite{Steel2008} a naive estimation of the $\lambda_{i\cdot}$ scale parameters will lead to very
permissive values, and thus we employ an informative hierarchical exponential model for them -- that is, we explicitly
inform the model that we prefer small distances overall.
The $\lambda_{0j}$ variables themselves come from exponential distributions with a fixed parameter.  
Since distinct gene families will usually have a different number of members, we must ensure that scale parameters 
$\lambda_{ij}$ can be comparable across gene families.
This is achieved through the scaling parameters $m_{ij}$, which are equivalent to a standardization of all distances
$d_j(\cdot,\cdot)$ to the interval $[0,1]$. 
The choice of a common prior parameter $\lambda_{0j}$ shared across gene families helps to avoid overparameterization.
But more importantly, it assumes that the underlying biological process reflected by the distance $j$ -- like the
duplication rate, for instance -- acts on the whole genome, while at the same time the $\lambda_{ij}$ parameters for the
observed distances are particular to each gene family $i$.

\subsection{Sampling from the posterior distribution}
Sampling directly from the above hierarchical model is very hard, so we employ a Markov chain Monte Carlo
sampling \citep{Liu2001} where variables are updated one by one using Generalized Multiple-try Metropolis (GMTM) updates
\citep{Liu2000, Pandolfi2010}.
The GMTM algorithm for proposing an update of current state $x$ to a new state $y$ can be described as
\begin{enumerate}
\item Given current state $x$, draw $k$\ \ samples $y_1,\hdots,y_k$ independently from the proposal distribution $p(\cdot\mid~x)$
\item Calculate the weights $\overline{w}_i(y_i;x)=\frac{w_i(y_i;x)}{\sum_{l=1}^k w_l(y_l;x)}$ where $w_i(y_i;x)$ can be any
        bounded, positive function 
\item Select $y=y_j\in{y_1,\hdots,y_k}$ according to $\overline{w}_i(y_i;x)\quad i=1,\hdots,k$ and set
        $W_y=\overline{w}_j(y;x)$ (assuming $j$ is the chosen index)
\item Similarly, draw $x_i^*\sim p(\cdot\mid~y)$ for $i=1,\hdots, j-1, j+1,\hdots, k$, and let $x_j^*=x$ (that is,
        using same index $j$ as before). Calculate $\overline{w}_i(x_i^*;y)=\frac{w_i(x_i^*;y)}{\sum_{l=1}^k w_l(x_l^*;y)}$ and
        set $W_x=\overline{w}_j(x;y)$
\item Accept $y_j$ as new state with probability
\begin{equation}
        \min\left(1,\frac{\pi(y_j)p(x\mid y_j)}{\pi(x)p(y_j\mid x)}\frac{W_x}{W_y}\right)
\end{equation}
\end{enumerate}
where $\pi(y)$ is the posterior distribution under $y$. 
A common choice for the weight $w(y;x)$ is the posterior distribution $\pi(y)$ (terms independent from $y$ will cancel
out).
This is how the vector $\bm\lambda_0$ is updated. 

When updating the gene trees, as mentioned before, we avoid direct calculation of the joint probability of alignment $D_i$ and model $\bm\theta_i$, 
$P(\bm{D},\bm{\theta}\mid \bm{G})$, assumed to be independent between gene families, and that can be expanded as the product
of the phylogenetic likelihood and the prior of the substitution model \citep{Yang1997,Larget1999}
 \begin{equation} 
P(\bm{D},\bm{\theta}\mid \bm{G}) = \prod_{i=1}^N P(D_i,{\bm \theta_i}\mid G_i) = \prod_{i=1}^N P(D_i\mid{\bm \theta_i}, G_i) P(\bm \theta_i)
\end{equation}
We furthermore assume that for each gene family $i$ we have a posterior distribution of trees $P(G_i\mid D_i, {\bm\theta_i})$ (for example estimated by
some state-of-the-art phylogenetic inference program) and use the posterior frequency of tree $\hat{G}$ as a resampling weight
$g_i(\hat{G})= P(G_i=\hat{G}\mid D_i, {\bm\theta_i})\propto P(D_i,{\bm \theta_i}\mid\hat{G})$ \citep{smith1992}.
This procedure is similar to the two-stage MCMC approach implemented in BUCKy \citep{Ane2007}, and in practice means that
by using the resampling weights as the proposal distribution $p(\cdot\mid\cdot)$ in the GMTM updates the terms $P(D_i,{\bm \theta_i}\mid G_i)$
will cancel out.
Under this approach the substitution parameters $\bm\theta_i$ and the branch lengths are nuisance parameters and are neglected from the
analysis.  At every iteration we propose a new species tree by applying one or several branch swaps on the current one, and we also have a rerooting proposal that does not change the topological information. 

\subsection{The normalization constant}\label{sec:normal}
When updating the species tree $S$ or the distance scale parameters $\lambda_{i\cdot}$, we must also update the
normalization constant -- also called the partition function \citep{Murray2012,Chung2013}. 
It is the denominator of our multivariate exponential distribution, and ensures that this distribution sums up to one over all possible gene topologies with the same leafset. 
It must be taken into account if we want to interpret the penalties as probabilities, and neglecting it can even change
the relative score of species trees under some distances \citep{Bryant2008}.
This function changes with the species tree and $\lambda_{i.}$, and computing it is impractical for gene trees with more than a few leaves.
With the GMTM algorithm the weight function is not restricted to probability measures from the model.
Hence we could postpone the partition function calculation to step (5) of the GMTM algorithm by employing the unscaled
distance distribution as the weight function.
However the acceptance probability will depend on the ratio between the partition functions
$Z_i(S^*,\lambda_{i\cdot}^*)/Z_i(S,\lambda_{i\cdot})$ which do not cancel out.

Our solution is to use an exchange algorithm, where the proposal of a parameter that affects the partition function 
is always accompanied by the proposal of an auxiliary gene tree, therefore canceling out the partition functions
\citep{Murray2012,Liang2010,Atchade2008,Caimo2012}.
This auxiliary gene tree is an augmented variable which is unrelated to our resampled estimates -- it is sampled through
a secondary MCMC given the proposed $(S^*,\lambda_{i\cdot}^*)$ and is discarded immediately after calculating 
the exchange ratio \citep{Caimo2010,Caimo2012}.
        
Namely, given the unscaled distance distribution $f_i(G,S,\lambda_{i\cdot})=e^{-\sum_{j=1}^J d_j(G,S)/(m_{ij}\lambda_{ij})}$ 
then the joint distance probability can be written as  
\begin{equation}
        P(\bm{G}\mid \bm{\lambda},S) = 
        \frac{\prod_{j=1}^N f_i(G_i,S,\lambda_{i\cdot})}{\prod_{i=1}^N Z_i(S,\lambda_{i\cdot})} =
        \frac{f(\bm{G},S,\bm{\lambda})}{\prod_{i=1}^N Z_i(S,\lambda_{i\cdot})} 
\end{equation}
Then the GMTM exchange algorithm proceeds as following, 
remembering that each state is composed of the species tree $S$ and the vector $\bm{\lambda}$, but in practice we update
just $S$ or just a block of $\lambda_{i\cdot}$ for a gene family $i$ at a time. 
\begin{enumerate}
        \item Given current state $x$, draw $k$\ \ samples $y_1,\hdots,y_k$ independently from the proposal distribution $p(\cdot\mid~x)$
\item Calculate the weights $\overline{w}_i(y_i;x)=\frac{w_i(y_i;x)}{\sum_{l=1}^k w_l(y_l;x)}$ where $w_i(y_i;x)=f(\bm{G},S^{(i)},\bm{\lambda}^{(i)})$
        is the joint unscaled distribution 
\item Select $y=y_j\in{y_1,\hdots,y_k}$ according to $\overline{w}_i(y_i;x)\quad i=1,\hdots,k$ and set
        $W_y=\overline{w}_j(y;x)$
\item Draw $x_i^*\sim p(\cdot\mid~y)$ for $i=1,\hdots, j-1, j+1,\hdots, k$, and let $x_j^*=x$
        Calculate $\overline{w}_i(x_i^*;y)=\frac{w_i(x_i^*;y)}{\sum_{l=1}^k w_l(x_l^*;y)}$ and
        set $W_x=\overline{w}_j(x;y)$
\item Assuming that the proposed state is $y_j=(S^*,\bm{\lambda}^*)$, then draw a set of auxiliary gene trees
        $\bm{A}=(A_1,\hdots,A_N)$ from $P(\cdot\mid\bm{\lambda}^*,S^*)$. This can be achieved by MCMC starting
        at the current $\bm{G}$ (notice that within this MCMC the partition functions cancel out.)
\item Accept $y_j=(S^*,\bm{\lambda}^*)$ as new state with probability
\begin{equation}
 \min\left(1,\frac{P(\bm{\lambda}^*\mid\bm{\lambda_0})p(x\mid y_j)}{P(\bm{\lambda}\mid\bm{\lambda_0})p(y_j\mid x)}
 \frac{f(\bm{G}\mid\bm{\lambda}^*,S^*)f(\bm{A}\mid\bm{\lambda},S)}{f(\bm{A}\mid\bm{\lambda}^*,S^*)f(\bm{G}\mid\bm{\lambda},S)}
 \frac{W_x}{W_y}\right)
\end{equation}
\end{enumerate}

Our preliminary results suggested  that point estimates from the posterior distribution of species trees are not affected by
neglecting the partition function, at least when using multiple distances as we did. 
Therefore, to reduce the computational workload here we will show results where this normalization constant has been neglected.

\subsection{Maximum likelihood estimation of species trees}
We furthermore implemented a maximum likelihood version of our model, where instead of sampling from the posterior
distribution $P(S,\bm{\Theta}\mid\bm{D})$ we try to find only the optimal values $\hat{S}$ and
$\hat{\bm{\Theta}}$ by incrementally decreasing the temperature $kT$ of the distribution $P(S,\bm{\Theta}\mid\bm{D})^{1/kT}$.
Therefore we can use the same algorithms as before, since as the temperature $kT$ decreases only updates that improve
the posterior distribution are accepted \citep{Rubenthaler09}.
The final state of this modified MCMC sampler is then its simulated annealing estimate, which can be interpreted as a ML
supertree estimate under an arbitrary set of distances (of which DLI and DLIR have been employed in the present study,
as we will see next.)

\subsection{Software implementation: {\em guenomu}}
We implemented our Bayesian supertree model into the program {\em guenomu} (\url{http://bitbucket.org/leomrtns/guenomu}), where the user can 
choose which distances will be taken into account by the multivariate exponential
distribution. Below we work with two parameterizations: one considering only the reconciliation distances, which we will
call the DLI model since it employs the minimum number of Duplications, Losses and Incomplete lineage sortings; and another adding the $mulRF$ distance, called the DLIR model. 
The program {\em guenomu}  is paralellized at the gene family level, whereby the communication is minimized through the
implementation of shared pseudo-random generation streams \citep{Feng2003,Feng2006}.
It also allows for implementing other distances without compromising the model.

\subsection{Assesment of accuracy using computer simulations}
In order to assess the performance of our implementation, we carried out a simulation study to compare  {\em guenomu's}
posterior species tree distribution with the true species trees and with other similar supertree methods described
below. We generated 12348 replicate data sets using {\em SimPhy} \citep{Mallo2013}. Our simulation scenarios takes
into account the birth and death of (paralogous) loci inside a species tree as well as the coalescent process describing
the gene tree within each locus.

For each replicate data set,  a species tree was generated under the Kingman coalescent process using the software
Dendropy \citep{Sukumaran2010} such that the number of species was between $6$ and $100$ and the length from the root to
tips was between $10^2$ and $10^4$ generations. We assumed that the effective population size was constant  and that the
species tree length in coalescent units \citep{Leache2011} was between $0.05$ and $5$. Then we used {\em SimPhy}  to
simulate the evolution of a number of gene families inside each species tree, generating gene trees with branch lengths
in substitution units.  Input parameters for each replicate were sampled stochastically from the distributions
summarized in table \ref{tab:params}.
In short, for each gene family the number of individuals per species varied between $1$ and $10$ while the number of
gene families ranged from between $3$ and $50$. Furthermore each gene family differed from the species tree through a
birth-and-death process of duplications and losses such that the average number of duplications was between $10^{-3}$
and $4$. In the end, each replicate had gene trees with average sizes ranging from $19$ to $1430$ leaves. Both the
generated species tree and each gene family tree are ultrametric, with lenghts given in generation times. But we used
{\em SimPhy} to further simulate substitution rate heterogeneities at several levels, by applying gamma-distributed rate
multipliers at the branches of each of the simulated trees from the three-tree model, namely the species tree, the locus
tree and the gene tree \citep{Mallo2013} . As a result the gene family trees will not follow a molecular clock.  Although
our model for species tree inference do not take into account this branch length information, it is relevant to simulate
gene tree uncertainty (i.e., tree inference error), as we will see below.

\medskip 
[Table \ref{tab:params} around here] 
\medskip

\subsection{Gene tree uncertainty}
The simulation procedure just described  provide us with the single, true tree underlying each gene family. However, our
model explicity considers gene tree uncertainty, relying as input on a sample of the posterior distribution of trees for
each gene family. To obtain these we could have simulated sequences along the true gene trees
\citep{Yang2007,Fletcher2009}, and then infer posterior gene tree distributions using  popular Bayesian phylogenetics
programs like MrBayes \citep{Ronquist2012}, PhyloBayes \citep{Lartillot2004} or BEAST \citep{Drummond2012}. 
Unfortunately this approach would be extremely slow, limiting the scope of our simulations, so we devised a way of
directly building gene tree distributions for each gene family in an attempt to recreate gene tree uncertainty. 

For each true gene tree we generated a population of topologies that may differ from the original by one or more
nearest-neighbor interchanges \citep{Waterman1978}, assuming a probability $p_T$ of applying branch swaps to each tree,
and a probability of branch swap per branch inversely proportional to its length, bounded by the maximum swap
probability $p_B$. The latter accounts for the observation that shorter branches are harder to reconstruct, and
generally are where most of the uncertainty in the phylogenetic estimation is concentrated. Branch swaps implied at
random one of the two nearest neighbors. In this way, a proportion $(1-p_T)$ of resulting trees are expected to be
identical to the original gene tree, while the others should differ by at least one branch. Values for $p_B$ larger than
1 are allowed, meaning that even larger branches will have maximum uncertainty. For low levels of uncertainty the most
frequent topologies will match the true gene family trees, although for the parameter values employed in our simulations
almost all data sets contained at least one most frequent gene tree distinct from the true one.
The pseudo-code for applying uncertainty to a gene tree is described in the Appendix.

\medskip 
[ Figure \ref{fig:simul} around here ]
\medskip 

The general workflow of the simulation study  is represented in Figure  \ref{fig:simul}.
Each replicate data set was composed by a species tree and a set of 3 -- 50 gene
families, where  phylogenetic uncertainty around each gene family is represented by a tree collection with up to 160
distinct topologies. As a result of the introduction of uncertainty, in 8\% of simulations all gene tree topology
distributions showed as  most frequent tree the original, true tree. In the remaining 92\%  at least one distribution
had a gene tree more frequent than the true tree. From these, wrong trees had on average a cummulative frequency 17\%
higher than the true tree, although in all cases the true topology could be found amongst the input gene trees,
even if at a low frequency. These gene tree frequencies were assumed by {\em guenomu} to be the resampling
weights (posterior frequencies from independent analysis) of each tree.

\subsection{Comparison with other species tree approaches}
For each simulated data set two Bayesian {\em guenomu} runs were conducted, under both the DLI and DLIR parameterizations.
Each run generates a (posterior) distribution of rooted species trees, and their equivalent posterior
distribution of unrooted trees. 
We used as point estimates of this distribution both the {\it maximum a posteriori (MAP) tree} and the 50\% majority
rule {\it consensus tree} -- remembering that they contain only the topological information and no branch lengths. 
Furthermore, under the DLIR model we conducted another run of our model to estimate a {\it maximum likelihood} species
tree  using simulated annealing, as described above.

For the same simulated data sets, and for comparison purposes, we also employed other supertree-related models to reconstruct
the species tree, namely  Gene Tree Parsimony \citep{Bansal13} and coalescent distance matrix methods \citep{Liu2009},
also called coalescent summary statistics methods. Because in this study  we consider gene trees evolving inside a
specie tree under GDL and ILS, we did not consider species tree methods that can not deal with multi-labelled gene
trees, like the supermatrix approaches or consensus methods like those implemented in BUCKy \citep{Ane2007}.

Gene Tree Parsimony aims to find species tree that minimizes its total reconciliation cost  with a given collection of
gene tree topologies. Here we used the software {\em iGTP} \citep{Chaudhary2010}, which assumes that the disagreement
between gene trees and the species tree is either  due to duplications and losses, or  to incomplete lineage sorting. An
independent analysis was conducted under each assumption. {\em iGTP} allows for each gene tree to have a weight
associated with it, and therefore we  used the topology frequencies resulting after the introduction of phylogenetic
uncertainty (see above). However, for computational feasibility we limited the input gene trees to the most frequent
quartile.

We  implemented in-house several summary statistics, coalescent species tree reconstruction methods based on  gene-wise
distance matrices, as described in \cite{Helmkamp2012}. 
This class of methods start by creating a matrix of distances
between species for each gene phylogeny. Such matrices are then merged into a single one, from which a species tree is
built using a clustering algorithm. 
The four procedures implemented are  GLASS \citep{Mossel2008},  STEAC
\citep{Liu2009e},  SD \citep{Maddison2006}, and  MAC \citep{Helmkamp2012} methods, which differ on how they handle
contradicting distances for the same pair of species. 
This distinction stems from the fact that when there is more than one individual from the same species in a given gene
family, its tree will have more than one path for the same pair of species. 
In such cases the GLASS and
SD methods will choose the minimum distance between the pair as the element gene-wise distance matrix, while the MAC and
STEAC methods will take the average between all distances for the same pair of species. However when there are several
gene families each matrix might have a distinct distance between a given pair of species. In these case, the GLASS and
the MAC methods use the minimum value, while the SD and STEAC use the mean distance accross matrices. The species tree
is then reconstructed based on this final distance matrix, where SD and STEAC employ an UPGMA, and GLASS and MAC use a
single-linkage clustering algorithm \citep{Helmkamp2012}.

These distance matrix methods are grounded on the divergence times between species, and thus rely heavily on the
individual branch lengths. Moreover, these methods assume a single gene tree for each gene family. However we assume
that, as with the previous approaches, there is gene tree uncertainty. To make a fair comparison with {\em guenomu} and
iGTP we generated gene tree distributions as before but only used the most frequent gene topology, assuming that all
branches had a length of $1.0$. The latter corresponds to using the path length (number of nodes) between leaves as the
pairwise distances \citep{Steel1993}. In addition, and for comparison purposes, we also used the true gene trees with the
true branch lengths. 
From each of these methods we obtain a single estimate of the species tree.

We also did a preliminary analysis with the mulRF supertree model \citep{Chaudhary2013} as implemented in-house, but
pilot experiments showed very bad performance and we did not explored it further. 
Due to their computational complexity, we did not evaluate sequence-based approaches like *BEAST \citep{Heled2010} as
each simulated gene family had hundreds of taxa, on average. 
Indeed, it would have been very interesting to compare our model with the recently developed software PhylDog, which
implements a probabilistic model to estimate the maximum likelihood species tree given a set gene families \citep{Boussau2013}.  
Unfortunately this software
is very difficult to install and run, and we were not successful in designing a pipeline that could systematically
execute PhylDog against our simulated data sets.

\subsection{Performance measures}
%In our simulations we have the true rooted species tree and their equivalent true unrooted gene trees, all with branch
%lengths. 

We called {\it tree accuracy} to the number of splits (bipartitions or branches) from the true unrooted species tree that
were successfully recovered in the inferred species tree, divided by the total number of splits on the true unrooted
species tree.  This measure will go from zero when no true splits are recovered to one when the true
species tree topology is completely recovered. We also recorded the proportion of simulations where the true unrooted
species tree was perfectly recovered (i.e.,  the proportion of replicates where tree accuracy was one, or {\it true tree
recovery}.) In addition, from the analyses were we had a sample from the posterior distribution we  also calculated
the proportion of credible sets (at the 95\% level) that contained the true unrooted species, or {\it true tree coverage}.

\subsection{Empirical benchmark}

To test our model on real data, we downloaded all maximum likelihood gene family trees from the TreeFam database
\citep{Schreiber2013} and then pruned all data outside the 12 {\it Drosophila} species subtree. We considered all the
gene families with more than 3 species.  Within these we identified 4591 single-copy gene families (that is, just
orthologs), and from the remaining 2562 gene families with at least one paralog, we arbitrarily selected the 43 largest
ones, with number of leaves between 102 and 295. As in the simulations, we added phylogenetic uncertainty around the
TreeFam gene families as decribed above, generating a distribution of gene trees for each gene family.

\section{Results}
In the simulation scenarios considered, our Bayesian supertree implementation performed best regardless of its
parameterization or tree point estimate, showing median tree accuracies around 0.7  
(Fig. \ref{fig:compareall}). 
Our Bayesian method was however followed very closely by iGTP, which was slightly more accurate when DL
costs instead of ILS costs were used.
Our ML estimates were almost identical to iGTP.
The GLASS method ranked next when given the true gene trees, but performed badly
when gene trees were uncertain.  
SD methods did worse, performing equally well on true and observed gene family trees.
By far, the worst methods were STEAC and in particular MAC, even when given the true gene trees.

\medskip 
[ Figure \ref{fig:compareall} around here ]
\medskip 

To better understand the effect of different parameter values used in the simulation we fitted a multiple linear
regression model to the tree accuracies obtained with {\em guenomu}'s consensus tree under the DLIR model
(Fig. \ref{fig:corr}). 
Several  explanatory variables were significantly and positively correlated with tree accuracy  
(p-value $<10^{-15}$ unless otherwise stated), including the smallest branch length and the sum of species tree branch lengths
in coalescent units, the expected number of duplications per gene family on the species tree, the average number of
individuals from the same species per gene family, the average number of distinct species per gene family (p-value
$<10^{-5}$) and the average gene family size (p-value $<10^{-6}$) . The species tree height (length from root to tips)
was also significantly correlated with tree accuracy according to simple linear regression, but this relationship faded
away once we included other variables like the smallest
branch and the total sum of branch lenghts. In any case, the most influential parameter was apparently the number of
species, which was negatively correlated with tree accuracy (p-value$<10^{-15}$).

% result from linear model  lm((1-a$V46) ~ ...)
%#  7 ... 9  : sum_sptreelengths, expected number of dups, observed number of species
%# 10 ... 12 : avge number indiv/species, avge n_species, avge n_leaves for gene tree [[avge over locus trees]] 
%# 13 ... 15 : gene tree inference error simulation (noise): avge prob of wrong trees (higher than true), error rate per tree and per branch, acc to bmc2_addTreeNoise
%a$V1         2.281e-03  7.293e-05  31.271  < 2e-16 ***
%tne          1.321e+00  8.153e-02  16.198  < 2e-16 ***
%TNe         -5.993e-03  5.207e-03  -1.151    0.250    
%sbe          1.553e-02  1.132e-03  13.724  < 2e-16 ***
%a$V8         1.498e-02  1.266e-03  11.830  < 2e-16 ***
%a$V9        -5.351e-03  2.347e-04 -22.794  < 2e-16 ***
%a$V10        4.106e-02  1.991e-03  20.620  < 2e-16 ***
%a$V11        1.249e-03  2.555e-04   4.891 1.02e-06 ***
%a$V12        6.947e-05  1.334e-05   5.207 1.95e-07 ***
%a$V13       -9.636e-03  1.241e-02  -0.777    0.437    
%a$V14       -5.806e-02  3.860e-02  -1.504    0.133    
%a$V15        1.112e-02  1.087e-02   1.023    0.306    

\medskip 
[ Figure \ref{fig:corr} around here ]
\medskip 

Neither {\em guenomu} nor {\em iGTP} were able to recover the true species tree for more
than 35 species (Fig. \ref{fig:true}).
Under guenomu, the true tree could be found within the 95\% coverage on data sets with up to 40 species.
Nonetheless, even when the true species tree could not be found, {\em guenomu} found trees consistently
more accurate than {\em iGTP} for large trees (Fig. \ref{fig:recover}). 
In these cases the consensus tree from the posterior samples under {\em guenomu} was in any case the best estimate. 
For example, for the simulations with more than 60 species, 58\% of the species trees reconstructed from {\em guenomu}
found at least 60\% of the true splits, against only 51\% of the trees inferred by {\em iGTP} (Fig. \ref{fig:recover}). 
5.7\% of trees from guenomu reconstructed more than 80\% of the true splits for the same large data sets, against 3.4\%
when {\em iGTP} was used.

\medskip 
[ Figure \ref{fig:true} around here ]
\medskip 

Finally, we fit another multiple linear regression model over all replicates where the dependent variable was the
difference in tree accuracy between {\em guenomu}'s consensus estimate from the DLIR model and {iGTP} under the DL cost.
This  indicated that {\em guenomu} significantly outperformed {\em iGTP} for higher levels of gene tree uncertainty as
measured by $p_T$ and $p_B$, as well as for lower duplication rates and more and larger gene families (at a p-value
lower than $10^{-7}$ each). The improvement of {\em guenomu} over {\em iGTP} was only marginally significant (p-value $<0.05$)
for larger species trees with shorter branch lenghts.

\medskip 
[ Figure \ref{fig:recover} around here ]
\medskip 

Besides the posterior distribution of species trees, {\em guenomu} also obtains for each replicate the posterior
distribution of gene trees, which correspond to reweighted versions of the input gene tree distributions. 
Gene tree posterior distributions were always more compact than the input gene tree distributions, with fewer trees having a
larger contribution, and this reduction was more pronounced for the DLIR parameterization than for DLI (results not
shown).

\subsection{TreeFam Drosophila gene families} 
As for the TreeFam {\it Drosophila} data, {\em guenomu} reconstructed the known species tree \citep{Stark2007} as well as
the root location for the data set with few large gene families (43 genes with more than 100 members each). The root
location could be found even in the absence of an outgroup.
The rooted species tree estimation from this data set was robust to gene tree uncertainty.
That is, even after replacing each ML tree by a distribution of topologies as described before, {\em guenomu} was capable
of finding the well-known species tree with posterior probability equal to one.

However, for the large data set composed of single-copy gene families (4591 gene families with up to 12 members each),
although the unrooted topology could be found in the absence of gene uncertainty, the rooted location was not correctly
inferred. 
In the presence of gene uncertainty, then even the reconstructed topology missed one or two splits.

\section{Discussion} 
\subsection{Reconstructing species trees from dynamic gene families with lineage sorting} 
Most phylogenomic analyses start with the identification of ortholog sequences, after some process of sequence
collection and pruning of paralogs (e.g. \cite{Koonin2008, Dunn2008, Medina2011, Sanderson2011,
Bininda2011,Williams2012,Lang2013, Romiguier2013, Salichos13}.)
Much less often, alternative species tree analyses focus on gene duplication and loss \citep{Katz2012}, or consider
paralogs at all \citep{Holton2010}.
The common practice of removing data (sequences or whole genes) before analysing the reduced data sets is sometimes
needed, in order to conform to model assumptions or due to computational limitations of the methods.
This process of data carving is also rooted on the elusive idea of removing uncertainty and improving resolution, which is
anathema to the Bayesian paradigm. Here we propose a Bayesian supertree model that does not need to separate
homolog from paralogs, that can work with multiple gene families and that is also able to deal with multiple individuals per
species.
It also minimizes the pre-processing of data to conform to strict standards of resolution, since the model can cope with
uncertainty in the phylogenetic inference.
Therefore our approach, implemented in the program {\em guenomu}, is able to work with all the available
data at hand. 
 
Our simulations suggest that {\em guenomu} is able to reconstruct reasonably accurate species trees
in the presence of gene duplications, losses and deep coalescences, even though it uses simple distances that
do not take into account branch lengths.
Although it is not a panacea, it suggests that exploring these apparently na\"ive models of disagreement may be an
attractive alternative to biologically realistic frameworks, specially once we aim for a consilience of methods.
We hope that more analyses are conducted under broadly usable frameworks like the one we offer, at least as a starting
point before embarking into other ones with additional assumptions. 

Both DLI and DLIR parameterizations performed almost identically in recovering the species tree, but the DLIR
model, which includes the mulRF distance \citep{Chaudhary2013}, helped reducing the uncertainty from the posterior gene
tree distributions. 
Since the computational overhead of the DLIR parameterization over DLI is minimal, it seems worthy to use DLIR.

In the  simulations the true species tree became more and more difficult to recover as the number of species increased.
However, even in these cases the majority of true splits were recovered through the consensus or MAP estimates,
suggesting that our model can be applied to relatively large species trees. We notice however that for the largest species
trees the true species tree was not covered in the 95\% credible interval, a behaviour expected from most phylogenetic
methods.

Surprisingly, iGTP also provided good species tree estimates, being only marginally worse than {\em guenomu} as the number of 
species increased. 
We did not expect this behaviour because iGTP cannot handle both DL and ILS costs at once, and nonetheless assuming just a DL
cost seemed to work quite well even in the presence of ILS.
Indeed, iGTP is  faster than {\em guenomu}, which then makes it a very competitive approach.
On the other hand, {\em guenomu} can replicate completely {\em iGTP}'s functionality through its ML algorithm, which could
furthermore work with other complex distance costs in the future.

The distance matrix methods evaluated here assume that the source of disagreement comes from the coalescent process,
therefore it comes at no surprise that they did not perform well under our simulations, with an important contribution
from duplications and losses. Duplications and losses, when unnacounted for, will have a deleterious effect on the
estimated coalescent process \citep{Rasmussen2012}. Although this comparison might seem somehow unfair to these summary
statistics methods, it serves also as a warning against using attractive algorithms without considering its assumptions.
It might be interesting, however, to explore more extensively these methods and see how they fare against {\em guenomu}
or {\em iGTP} when there is only ILS.

We also observed that both the GTP and our multivariate distance models benefit from large gene families, with several
representatives from each species -- even though we do not distinguish paralogs or populational samples. This testifies
against the trend of finding good, hand-crafted sets of orthologs at the expense of the paralogs that must be discarded
from the analysis. This pruning may be done after the complete data set was used to create the posterior distribution, in a subsequent
analysis, but not before.

%\citep{DeGiorgio2013}.
% We observed that the hyperprior $\lambda_{0\cdot}$ can influence the posterior distribution, with lower values
% producing  sharper distributions of species trees, which indicates that the data itself usually does not have enough
% information about the penalties under this model. Therefore we must always choose values small enough, informing the
% model about the allowed gene/species trees discrepancy.
%red{Here we have evaluated xx methods. In this comparison}

As expected, the length in coalescent units of the smallest branch of the species tree, and the overall species tree
length had a great influence on species tree accuracy according to our model and also for GTP methods, together with the
number of species.
Indeed, we expect more incomplete lineage sorting and consequently more disagreement between gene and species trees next
to short branches, and also over a short generation span within large effective population sizes.
Under the conditions simulated, our model performed better on data sets with more and larger gene families, but representing fewer species. 
Small gene families over a large number of species proved more challenging to the methods compared. 

The root location of the {\it Drosophila} largest gene families was found even without an outgroup since the presence of
paralogs is informative about the optimal location of duplication nodes \citep{Katz2012}.
For the single-copy genes data set, on the other hand, the root location was not correctly
inferred for the same reason: orthologous data sets lack such phylogenetic information.
Furthermore the smaller gene families were more sensitive to the effect of tree uncertainty.

\subsection{Advantages and shortcomings of our Bayesian supertree method}
Bayesian methods in general give us information about the reliability of their findings through the posterior distribution. 
The flatness of the posterior distribution of species trees, or the frequency of the most frequent topology can give us
a hint about the amount of signal contained in the data.
The use of a Bayesian framework not only allows for the test of hypothesis and for model selection, but also facilitates the
creation of models of increased complexity, subject always to computational constraints.
Thus our framework can be extended in a straightforward manner to handle more complex situations, perhaps incorporating
more biological realism. Our distinct parameterizations allow us to compare models that use different sets of distances,
such that we can, in principle, compare their underlying assumptions.

Our modelling of "phylogenomic error" can be seen as a generalization of optimal Bayes estimators, like the consensus
tree or quartet puzzling \citep{Huggins2011}.
Such a model has already been suggested as an alternative to the detailed description of the processes causing the
tree incongruences, while at the same time maintaining the convenience of supertree methods \citep{Steel2008, Cotton2009}.
It can help bridge the gap between sequence-based alignment and phylogenomic analyses \citep{Cotton2009}, and ours is
the first implementation of this class of methods.

Originally, parsimonious reconciliation methods under the DL or the ILS costs worked with rooted gene trees and rooted
species trees \citep{Guigo1996,Than2009}. 
However, since we define the D, L and ILS distances as the minimum cost over all possible gene root locations, our model
works with unrooted input gene trees.
This is an important advantage, as in general there is no reliable information about the rooting at the gene tree level, and fixing a
root location can severely influence the distances.

The current implementation of our model uses only the topological information, ignoring branch lengths, but it can be
extended in the future to deal explicitly with branch lengths on the gene and species trees. However, it is possible to
use crudely the branch lengths as an indication of confidence in the bipartition, as implemented in our generator of
gene tree uncertainty.

Other topological distances can be easily implemented as well, taking advantage of our solution to the normalization
constant.
An example of a promising distance would be one based on the recently developed DLCpar algorithm \citep{Wu2013}, which
can find the most parsimonious reconciliation scenario by considering duplications, losses, and deep coalescences at
once.
Another option would be based on algorithms for finding the minimum cost reconciliation in the presence of duplications,
losses and HGT \citep{Tofigh2011, Doyon2011}.
One drawback is that many of these algorithms do not return a single solution, since under these scenarios the best
reconciliations are not unique.

When using resampling weights it is important to have a large sample with all reasonable trees represented for each gene
family \citep{smith1992}, specially since even lower likelihood trees might contribute more to the species tree than the
ML estimate \citep{Boussau2013}. 
Therefore ideally under our importance sampling algorithm we should avoid using point estimates of the gene trees (like
the ML gene tree, for instance), although for the TreeFam data set it did not hamper the analysis.
Nonetheless this simplification allows us to focus on the distribution of tree distances while using well-established
software to estimate the individual gene tree distributions, saving time and effort.

\subsection{Acknowledgments}
We want to thank Nicolas Lartillot for his suggestions on the solution for the partition function. This research was
financially supported by the European Research Council (ERC- 2007-Stg 203161-PHYGENOM). DM was also supported by an FPI
fellowship (BES-2010-031014) from the Spanish Government at the University of Vigo.


%bibliography
\begin{spacing}{1.4} \bibliography{bibtex_probpars} \end{spacing}
%\begin{spacing}{0.8} \printbibliography \end{spacing}
\appendix
\section*{Appendix}
\renewcommand{\thesubsection}{\Alph{subsection}}
\subsection{pseudo-code for the generation of tree uncertainty}
Here are the steps to apply uncertainty to a gene tree $T$:
\small
\begin{enumerate}
\item find maximum branch length $t_{MAX}=MAX_i(t_i)$ over all internal branches $i$ of $T$, where $t_i$ is the branch length; 
\item associate to each internal branch $i$ a probability of swap $p_i=p_B\times (1 - t_i/t_{MAX})$;
\item with probability $(1-p_T)$, create a copy of $T$ unchanged; 
\item Otherwise, with probability $p_T$, apply uncertainty as follows. For each internal branch $i$ in postorder (that is, closer to leaves first), using arbitrary root location and assuming left and right children are labelled $c1$ and $c2$:
\begin{enumerate}
\item If $c1$ is internal but $c2$ is not, then 
\begin{itemize} 
\item with probability $p_{c1}/3$ attach the left child of $c1$ to $i$; 
\item otherwise with probability $p_{c1}/3$ attach the right child of $c1$ to $i$; 
\item otherwise do nothing;
\end{itemize}
\item If $c2$ is internal but $c1$ is not, then 
\begin{itemize} 
\item with probability $p_{c2}/3$ attach the left child of $c2$ to $i$; 
\item otherwise with probability $p_{c2}/3$ attach the right child of $c1$ to $i$; 
\item otherwise do nothing;
\end{itemize}
\item If both $c1$ and $c2$ are internal branches, then do one of the following:
\begin{itemize}
\item with probability $p_{c1}(1-p_{c2})$ attach the left or right child of $c1$ to $i$, at random;
\item with probability $p_{c2}(1-p_{c1})$ attach the left or right child of $c2$ to $i$, at random;
\item with probability $p_{c1}p_{c2}/2$  attach the left or right child of $c1$ to $i$, then attach the left or right child of $c2$ to one of $i$, the left or right children of $c1$;
\item with probability $p_{c1}p_{c2}/2$  attach the left or right child of $c2$ to $i$, then attach the left or right child of $c1$ to one of $i$, the left or right children of $c2$;
\item with probability $(1-p_{c1})(1-p_{c2})$, do nothing;
\end{itemize}
\end{enumerate}
\end{enumerate}
\normalsize
If the gene tree does not have branch lengths, then a common value is assumed for all branches.
Furthermore if there are several trees in a file, the above procedure is repeated for each one of them, keeping track of
their frequencies.

\newpage
\noindent{\large \bf Tables} 
\begin{table}[htbp]
\begin{spacing}{1.4} 
\caption{\label{tab:params} Parameter values used in the simulations. Each simulation replicate was parameterized with values sampled from predefined statistical distributions. Given a species tree, simulated with {\em Dendropy}, the program {\em SimPhy} was used to simulate  gene family trees. An inhouse program was then used to mimic gene tree uncertainty, transforming each gene tree in a distribution of
gene tree topologies in an attempt to emulate the effect of alignment-based phylogenetic inference in practice.
}
\begin{tabular}{c l c l}
        & \bf{description} & \bf{symbol} & \bf{distribution} \\
\cline{2-4} 
\bf{species tree} & & & \\
{\em (Dendropy)} & Number of species & & Uniform(6,$10^2$) \\
                  & Number of generations (total tree height) &  & Uniform($10^2$,$10^4$)\\
                  & Expected number of duplications & $E_{dup}$ & Uniform($10^{-3}$, 4)\\
                  & Number of gene families & $N$ & Uniform(3,50)\\
                    & Rate heterogeneity multiplier & $H_s$ & Gamma(1,1)\\
\cline{2-4} 
\bf{locus tree}${}^a$ & & & \\
{\em (SimPhy)}& Gene duplication rate ${}^b$ & $\beta$ & Exponential($E_{dup}/\sigma$)\\
                      & Gene loss rate & & Uniform(0, $0.75\times\beta$)\\
                        & Rate heterogeneity multiplier & $H_l$ & Gamma(1,1)\\
\cline{2-4} 
\bf{gene tree}${}^a$ & & & \\
{\em (SimPhy)}& Effective population size & $N_e$ & $2000 \times$ LogNormal(0,0.25)\\
                     & Number of individuals per species ${}^c$ & & Uniform(1,10)\\
                     & Substitution rate (per time unit) &&$0.001$\\
                     & Rate heterogeneity multiplier & $H_g$ & Gamma(1,1)\\
\cline{2-4} 
\bf{gene tree} & & & \\
\bf{uncertainty} & Maximum number of generated trees && 160\\
{(inhouse)} & Tree dispersion term & $L_T$ & Uniform(2,5)\\
                & Tree location term & $D_T$ & Uniform(3,6)\\
                & Frequency of trees with uncertainty & $p_T$ & $1.1 \times$ Beta($L_T\times D_T$,$D_T$)\\
                & Branch location term & $L_B$ & Uniform(1,5)\\
                & Branch-wise error probability & $p_B$ & $1.5 \times$ Beta($L_B$,1)\\ 
\hline
\multicolumn{4}{l}{\footnotesize ${}^a$ terminology used by {\em SimPhy}}\\
\multicolumn{4}{l}{\footnotesize ${}^b$ the term $\sigma$ is the total sum of branch lengths (in generations) in the species tree}\\
\multicolumn{4}{l}{\footnotesize ${}^c$ common for all gene families, emulating perfect sampling of individuals}\\
\end{tabular}
\end{spacing}
\end{table}
\newpage
\noindent{\large \bf Figure Captions} 

\begin{figure}[h]\centering
\begin{spacing}{1.66} 
%{\begin{minipage}{7.2in}\includegraphics[width=7.2in,height=3.in]{fig1.pdf}\end{minipage}}
\caption{\label{fig:simul}
Simulation workflow. On the left a single data set is produced using  {\em Dendropy} and {\em SimPhy}. The true species
tree (rooted, with branch lengths) simulated by {\em Dendropy} is used as input by {\em SimPhy} to generate several
(rooted) phylogenies, one per gene family. Then, uncertainty is added such that we have a distribution of topologies
(unrooted, no branch lengths) per gene family. These collections will be used as input to the inference programs in
the right panel. 
Our software {\em guenomu} estimates the posterior distribution of species trees (rooted, without
branch lengths) and the posterior distributions of gene family topologies, based on all input gene tree distributions.
The software {\em iGTP} also uses the input gene trees after transforming the frequencies into integer values
representing weights, and estimates two rooted species trees: one under the Duplication and Loss cost, and another under
the Deep Coalescence cost. For the distance-based species tree inference algorithms only one tree per gene family is
used, and two alternative choices were attempted: one was to use the true gene families, with branch lengths; and the
other was to use the most frequent 
gene trees (topologies only) after introducing phylogenetic uncertainty. 
In the later case it was assumed that all branches had the same length.
}\end{spacing}\end{figure}

\begin{figure}[h]\centering
\begin{spacing}{1.66} 
%{\begin{minipage}{7in}\includegraphics[width=7in]{fig2.pdf}\end{minipage}}
\caption{\label{fig:compareall}
Tree accuracy. 
Each violin plot (kernel density curves plus a boxplot) represents the distribution of
tree accuracy values for the different methods evaluated. 
We ran {\em guenomu} 3 times for each simulation, to sample from the posterior distribution under  the DLI and
the DLIR models ('DLI' and 'DLIR'), and  to find the ML estimate under the DLIR parameterization ('DLIR-ML'). 
From each posterior distribution we
obtained two point estimates of the species tree: the MAP tree and the consensus tree (labelled 'MAP' and
'cons'). 
We also obtained two estimates of species trees using iGTP, running the program under the DL cost and
under the ILS cost ('iGTP-DL' and 'iGTP-ILS'). We also ran 4 distance matrix approaches (GLASS, SD, STEAC and MAC) using
two types of data sets, just the topologies with uncertainty ('topo') or the true simulated gene family trees with
branch lengths ('true'). 
At the top we show the hierarchical clustering of the different methods based on their tree accuracy values for all
replicates.
}\end{spacing}\end{figure}

\begin{figure}[h]\centering
\begin{spacing}{1.66} 
%{\begin{minipage}{7.1in}\includegraphics[width=7in]{fig3.pdf}\end{minipage}}
\caption{\label{fig:corr}
{\em Guenomu}'s tree accuracy with respect to simulation parameters.
Here the species tree was estimated as the consensus tree from the posterior distribution assuming a DLIR
parameterization.
Based on a multiple linear regression analysis, these are the parameters
that most significantly affected tree accuracy: at the top we have, from left to right, the total number of
gene families, the total species tree length (in coalescent units), the total number of species on the species tree, and
the average gene family size (number of leaves); at the bottom we have the length of the smallest branch in the species
tree (in coalescent units), the expected number of duplications per replicate, the average number of individuals from
same species per gene family, and the average number of species represented by gene family. Over each panel we show the
regression line overlaid in red.
}\end{spacing}\end{figure}

\begin{figure}[h]\centering
\begin{spacing}{1.66} 
%{\begin{minipage}{5in}\includegraphics[width=5in]{fig4.pdf}\end{minipage}}
\caption{\label{fig:true}
Tree accuracy and true tree recovery for {\em guenomu} and {\em iGTP} according to the number of species.
The lower lines show the true tree recoveries pooled over the number of species, while the higher lines show the
average tree accuracy under different conditions.
The gray area represents the true tree coverage under {\em guenomu} DLIR model (ie., the proportion of data sets where the true species tree was 
within the 95\% credible set).
Each line is actually a smooth regression line over the individual values. 
For {\em iGTP} the species tree was estimated under the DL cost, while for {\em guenomu} we plot the ML estimate as well
as the Bayesian consensus, both under the DLIR parameterization.
}\end{spacing}\end{figure}

\begin{figure}[h]\centering
\begin{spacing}{1.66} 
%{\begin{minipage}{5in}\includegraphics[width=5in]{fig5.pdf}\end{minipage}}
\caption{\label{fig:recover}
Tree accuracies of {\em iGTP} and {\em guenomu} for large data sets.
On the left we have all data sets where the number of species was between 30 and 60, and on the right we have the pooled
replicates with species tree sizes between 61 and 100 leaves. 
Each bar represents the proportion of data sets with tree accuracies between the ranges described in the middle, such
that  at the top we have the fraction of simulations where more than 80\% of the splits were recovered, while the
inferences where less than 20\% of the true splits were found are the bottom.
}\end{spacing}\end{figure}


\end{document}
````

## Source: paper001_guenomu/bibtex_probpars.bib
SHA-256: `5782fb5b53ffa8b358843712d6bd278c9a9d46c3fefaeab03c19d127bbc81170`

````bibtex
@BOOK{gelman2003,
  title = {Bayesian Data Analysis, Second Edition},
  publisher = {Chapman \& Hall/CRC},
  year = {2003},
  author = {Andrew Gelman and John B. Carlin and Hal S. Stern and Donald B. Rubin},
  pages = {696},
  edition = {2},
  month = jul,
  isbn = {158488388X}
}

@book{olive2005,
  title = {Applied Robust Statistics},
  author = {David Olive},
  year = {2005},
  publisher = {preprint},
  url = {http://www.math.siu.edu/olive/ol-bookp.htm}
}

@article{al2008,
  title={{Estimation of the Mean of Truncated Exponential Distribution}},
  author={Al-Athari, F.M.},
  journal={Journal of Mathematics and Statistics},
  volume={4},
  number={4},
  pages={284--288},
  year={2008}
}

% Multiple-try Metropolis
@Article{Liu2000,
  Title                    = {{The multiple-try method and local optimization in Metropolis sampling}},
  Author                   = {Liu, JS and Liang, F and Wong, WH},
  Journal                  = {Journal of the American Statistical Association},
  Year                     = {2000},
  Pages                    = {121},
  number = {449},
  volume = {95},
  publisher = {American Statistical Association},

  File                     = {:home/leo/articles/02-ImportedToMendeley/JAmStatAssoc.2000.121.pdf:pdf},
  Keywords                 = {MCMC},
  Mendeley-tags            = {MCMC},
  Url                      = {http://www.tandfonline.com/doi/abs/10.1080/01621459.2000.10473908}
}
@article{Craiu2007,
author = {Craiu, Radu V. and Lemieux, Christiane},
doi = {10.1007/s11222-006-9009-4},
file = {:home/leo/articles/MCMC.method/gibbs\_metropolis/StatsComput.2007.109.pdf:pdf},
issn = {0960-3174},
journal = {Statistics and Computing},
keywords = {antithetic variates,carlo,extreme antithesis,korobov rule,latin hypercube,markov chain monte,multiple-try metropolis,quasi monte carlo,random-ray monte carlo,sampling,sequence,sobol},
number = {2},
pages = {109--120},
title = {{Acceleration of the Multiple-Try Metropolis algorithm using antithetic and stratified sampling}},
url = {http://www.springerlink.com/index/10.1007/s11222-006-9009-4},
volume = {17},
year = {2007}
}
@article{Craiu2009,
author = {Craiu, Radu V.},
doi = {10.1198/jasa.2009.tm08393},
file = {:home/leo/articles/MCMC.method/adaptive\_parallel.mcmc/JAmStatAssoc.2009.1454.pdf:pdf},
issn = {0162-1459},
journal = {Journal of the American Statistical Association},
keywords = {adaptive markov chain monte,carlo,metropolis sampling,parallel chains,random walk metropolis sampling,regional},
number = {488},
pages = {1454--1466},
title = {{Learn From Thy Neighbor: Parallel-Chain and Regional Adaptive MCMC}},
url = {http://pubs.amstat.org/doi/abs/10.1198/jasa.2009.tm08393},
volume = {104},
year = {2009}
}
@article{Eidsvik2006,
author = {Eidsvik, Jo and Tjelmeland, H\aa Kon},
doi = {10.1007/s11222-006-5536-2},
file = {:home/leo/articles/MCMC.method/rj.mcmc/StatsComp.2006.93.pdf:pdf},
issn = {0960-3174},
journal = {Statistics and Computing},
keywords = {adaptive direction sampling,angular gaussian distribution,auxiliary variables,carlo,hastings,markov chain monte,metropolis,reversible jump mcmc},
number = {1},
pages = {93--106},
title = {{On directional Metropolis–Hastings algorithms}},
url = {http://www.springerlink.com/index/10.1007/s11222-006-5536-2},
volume = {16},
year = {2006}
}
@article{Andrieu2008,
author = {Andrieu, Christophe and Thoms, J},
file = {:home/leo/articles/MCMC.method/adaptive\_parallel.mcmc/StatsComput.2008.343.pdf:pdf},
journal = {Statistics and Computing},
number = {4},
pages = {343–373},
publisher = {Springer},
title = {{A tutorial on adaptive MCMC}},
url = {http://www.springerlink.com/index/979087678366R78V.pdf},
volume = {18},
year = {2008}
}
@Article{Pandolfi2010,
  Title                    = {{A generalization of the Multiple-try Metropolis algorithm for Bayesian estimation and model selection}},
  Author                   = {Pandolfi, Silvia and Bartolucci, Francesco and Friel, PN},
  Journal                  = {AISTATS},
  Year                     = {2010},
  Pages                    = {581--588},
  Volume                   = {9},

  File                     = {:home/leo/articles/02-ImportedToMendeley/AISTATS.2010.581.pdf:pdf},
  Keywords                 = {MCMC},
  Mendeley-tags            = {MCMC},
  Url                      = {http://eprints.pascal-network.org/archive/00007723/}
}
% MTM with different proposals
@article{Casarin2010,
abstract = {We propose a new class of interacting Markov chain Monte Carlo (MCMC) algorithms designed for increasing the efficiency of a modified multiple-try Metropolis (MTM) algorithm. The extension with respect to the existing MCMC literature is twofold. The sampler proposed extends the basic MTM algorithm by allowing different proposal distributions in the multiple-try generation step. We exploit the structure of the MTM algorithm with different proposal distributions to naturally introduce an interacting MTM mechanism (IMTM) that expands the class of population Monte Carlo methods. We show the validity of the algorithm and discuss the choice of the selection weights and of the different proposals. We provide numerical studies which show that the new algorithm can perform better than the basic MTM algorithm and that the interaction mechanism allows the IMTM to efficiently explore the state space.},
archivePrefix = {arXiv},
arxivId = {1011.1170},
author = {Casarin, Roberto and Craiu, Radu V. and Leisen, Fabrizio},
doi = {10.1007/s11222-011-9301-9},
eprint = {1011.1170},
file = {:home/leo/articles/02-ImportedToMendeley/arxiv.1011.1170v1.pdf:pdf},
journal = {arxiv},
month = nov,
number = {},
pages = {1--29},
title = {{Interacting Multiple Try Algorithms with Different Proposal Distributions}},
url = {http://arxiv.org/abs/1011.1170},
year = {2010}
}
%G-MTM
@article{Martino2012,
abstract = {The Multiple Try Metropolis method is an extension of the classical Metropolis algorithm in which the next state of the chain is chosen among a set of samples, according to normalized weights. In literature, a generalization of this technique has been proposed in Pandolfi et al. (2010) where the weight functions are not analytically specified. In this work, we propose different Multiple Try Metropolis schemes where the analytic form of weight functions can be chosen arbitrarily. The resulting algorithms fulfill the detailed balance condition.},
archivePrefix = {arXiv},
arxivId = {1201.0646},
author = {Martino, Luca and Read, Jesse},
eprint = {1201.0646},
file = {:home/leo/articles/02-ImportedToMendeley/arxiv.1201.0646v1.pdf:pdf},
isbn = {0034916249192},
journal = {arxiv},
keywords = {Bayesian,Computation,MCMC,MTM,Methodology,Statistics,Theory,algorithm},
mendeley-tags = {Bayesian,MCMC,MTM,algorithm},
month = jan,
title = {{Different acceptance functions for Multiple Try Metropolis schemes}},
url = {http://arxiv.org/abs/1201.0646},
year = {2012}
}

% normalization constant (doubly intractable models)
@Article{Atchade2008,
  Title                    = {Bayesian computation for statistical models with intractable normalizing constants},
  Author                   = {Yves Atchade and Nicolas Lartillot and Christian P. Robert},
  Journal                  = {arxiv},
  Year                     = {2008},

  Month                    = apr,

  Abstract                 = {This paper deals with some computational aspects in the Bayesian analysis of statistical models with intractable normalizing constants. In the presence of intractable normalizing constants in the likelihood function, traditional MCMC methods cannot be applied. We propose an approach to sample from such posterior distributions. The method can be thought as a Bayesian version of the MCMC-MLE approach of Geyer and Thompson (1992). To the best of our knowledge, this is the first general and asymptotically consistent Monte Carlo method for such problems. We illustrate the method with examples from image segmentation and social network modeling. We study as well the asymptotic behavior of the algorithm and obtain a strong law of large numbers for empirical averages.},
  Archiveprefix            = {arXiv},
  Arxivid                  = {0804.3152},
  volume = {0804.3152},
  Comments                 = {20 pages, 4 figures, submitted for publication},
  Eprint                   = {0804.3152},
  File                     = {:home/leo/articles/021/arxiv.0804.3152v2.pdf:PDF},
  Keywords                 = {00.read,Gibbs sampling,MCMC,adaptive mcmc,and phrases,bayesian inference,image segmentation,intractable constant,ising model,monte carlo methods,social network modeling},
  Oai2identifier           = {0804.3152},
  Owner                    = {leo},
  Timestamp                = {20130929.0003}
}
@article{Caimo2012,
abstract = {Exponential random graph models are a class of widely used exponential family models for social networks. The topological structure of an observed network is modeled by the relative prevalence of a set of local sub-graph configurations termed network statistics. One of the key tasks in the application of these models is which network statistics to include in the model. This can be thought of as statistical model selection problem. This is a very challenging problem---the posterior distribution for each model is often termed "doubly intractable" since computation of the likelihood is rarely available, but also, the evidence of the posterior is, as usual, also intractable. We present a fully Bayesian model selection method based on a Markov chain Monte Carlo algorithm of Caimo and Friel (2011) which estimates the posterior probability for each competing model as well as a possible approach for computing the model evidence.},
archivePrefix = {arXiv},
Volume = {1201.2337},
author = {Caimo, Alberto and Friel, Nial},
eprint = {1201.2337},
file = {:home/leo/articles/02-ImportedToMendeley/arxiv.1201.2337v2.pdf:pdf},
journal = {arxiv},
keywords = {Gibbs sampling,MCMC,intractable constant},
mendeley-tags = {Gibbs sampling,MCMC,intractable constant},
month = jan,
pages = {23},
title = {{Bayesian model selection for exponential random graph models}},
url = {http://arxiv.org/abs/1201.2337},
year = {2012}
}
@article{Liang2010,
author = {Liang, Faming},
doi = {10.1080/00949650902882162},
file = {:home/leo/articles/02-ImportedToMendeley/doublyIntractDist.10.1.1.182.3259.pdf:pdf;:home/leo/articles/02-ImportedToMendeley/JStatComputSimul.2010.1007.pdf:pdf},
issn = {0094-9655},
journal = {Journal of Statistical Computation and Simulation},
keywords = {00.read,Gibbs sampling,MCMC,autologistic model,autonormal model,auxiliary variable mcmc algorithm,change algorithm,ex-,intractable constant,metropolis-hastings algorithm},
mendeley-tags = {00.read,Gibbs sampling,MCMC,intractable constant},
month = sep,
number = {9},
pages = {1007--1022},
title = {{A double Metropolis-Hastings sampler for spatial models with intractable normalizing constants}},
url = {http://www.tandfonline.com/doi/abs/10.1080/00949650902882162},
volume = {80},
year = {2010}
}
@article{Friel2010,
abstract = {This paper proposes a new probabilistic classification algorithm using a Markov random field approach. The joint distribution of class labels is explicitly modelled using the distances between feature vectors. Intuitively, a class label should depend more on class labels which are closer in the feature space, than those which are further away. Our approach builds on previous work by Holmes and Adams (2002, 2003) and Cucala et al. (2008). Our work shares many of the advantages of these approaches in providing a probabilistic basis for the statistical inference. In comparison to previous work, we present a more efficient computational algorithm to overcome the intractability of the Markov random field model. The results of our algorithm are encouraging in comparison to the k-nearest neighbour algorithm.},
archivePrefix = {arXiv},
arxivId = {1004.3925},
author = {Friel, Nial and Pettitt, Anthony N.},
doi = {10.1007/s11222-010-9179-y},
eprint = {1004.3925},
file = {:home/leo/articles/02-ImportedToMendeley/arxiv.1004.3925.pdf:pdf},
journal = {arxiv},
keywords = {00.read,Gibbs sampling,MCMC,intractable constant},
mendeley-tags = {00.read,Gibbs sampling,MCMC,intractable constant},
month = apr,
pages = {12},
title = {{Classification using distance nearest neighbours}},
url = {http://arxiv.org/abs/1004.3925},
year = {2010}
}
@article{Caimo2010,
abstract = {Exponential random graph models are extremely difficult models to handle from a statistical viewpoint, since their normalising constant, which depends on model parameters, is available only in very trivial cases. We show how inference can be carried out in a Bayesian framework using a MCMC algorithm, which circumvents the need to calculate the normalising constants. We use a population MCMC approach which accelerates convergence and improves mixing of the Markov chain. This approach improves performance with respect to the Monte Carlo maximum likelihood method of Geyer and Thompson (1992).},
archivePrefix = {arXiv},
Volume = {1007.5192},
author = {Caimo, Alberto and Friel, Nial},
eprint = {1007.5192},
file = {:home/leo/articles/02-ImportedToMendeley/arxiv.1007.5192.pdf:pdf},
journal = {arxiv},
keywords = {00.read,Gibbs sampling,MCMC,intractable constant},
mendeley-tags = {00.read,Gibbs sampling,MCMC,intractable constant},
month = jul,
pages = {29},
title = {{Bayesian inference for exponential random graph models}},
url = {http://arxiv.org/abs/1007.5192},
year = {2010}
}

@ARTICLE{Muller2006,
  author = {M\o ller, J and Pettitt, AN and Reeves, R and Berthelsen, KK},
  title = {{An efficient Markov chain Monte Carlo method for distributions with
  intractable normalising constants}},
  journal = {Biometrika},
  year = {2006},
  pages = {451--458},
  file = {:home/leo/articles/02-ImportedToMendeley/Biometrika.2006.451.pdf:pdf;:home/leo/articles/02-ImportedToMendeley/Biometrika.2006.451\_expanded.pdf:pdf},
  keywords = {Gibbs sampling,MCMC,intractable constant},
  mendeley-tags = {Gibbs sampling,MCMC,intractable constant},
  url = {http://biomet.oxfordjournals.org/content/93/2/451.short}
}
@ARTICLE{Murray2012,
  author = {Murray, Iain and Ghahramani, Zoubin and MacKay, David},
  title = {{MCMC for doubly-intractable distributions}},
  journal = {arxiv},
  year = {2012},
  month = jun,
  abstract = {Markov Chain Monte Carlo (MCMC) algorithms are routinely used to draw
  samples from distributions with intractable normalization constants.
  However, standard MCMC algorithms do not apply to doubly-intractable
  distributions in which there are additional parameter-dependent normalization
  terms; for example, the posterior over parameters of an undirected
  graphical model. An ingenious auxiliary-variable scheme (Moeller
  et al., 2004) offers a solution: exact sampling (Propp and Wilson,
  1996) is used to sample from a Metropolis-Hastings proposal for which
  the acceptance probability is tractable. Unfortunately the acceptance
  probability of these expensive updates can be low. This paper provides
  a generalization of Moeller et al. (2004) and a new MCMC algorithm,
  which obtains better acceptance probabilities for the same amount
  of exact sampling, and removes the need to estimate model parameters
  before sampling begins.},
  archiveprefix = {arXiv},
  arxivid = {1206.6848},
  volume = {1206.6848},
  eprint = {1206.6848},
  file = {:home/leo/articles/02-ImportedToMendeley/arxiv.1206.6848.pdf:pdf},
  keywords = {Gibbs sampling,MCMC,intractable constant},
  mendeley-tags = {Gibbs sampling,MCMC,intractable constant},
  url = {http://arxiv.org/abs/1206.6848}
}


%% K2P distance 

@ARTICLE{Yang1999,
  author = {Yang, Z. and Yoder, Anne D.},
  title = {{Estimation of the Transition/Transversion Rate Bias and Species Sampling}},
  journal = {{Journal of Molecular Evolution}},
  year = {1999},
  volume = {48},
  pages = {274-283},
  number = {3},
  keywords = {, positive selection},
  publisher = {{Springer}},
}

%% STAR, STEAC, GLASS, Maximum trees

@article{Liu2010,
abstract = {We propose a model based approach to use multiple gene trees to estimate the species tree. The coalescent process requires that gene divergences occur earlier than species divergences when there is any polymorphism in the ancestral species. Under this scenario, speciation times are restricted to be smaller than the corresponding gene split times. The maximum tree (MT) is the tree with the largest possible speciation times in the space of species trees restricted by available gene trees. If all populations have the same population size, the MT is the maximum likelihood estimate of the species tree. It can be shown the MT is a consistent estimator of the species tree even when the MT is built upon the estimates of the true gene trees if the gene tree estimates are statistically consistent. The MT converges in probability to the true species tree at an exponential rate.},
author = {Liu, Liang and Yu, Lili and Pearl, Dennis K},
doi = {10.1007/s00285-009-0260-0},
file = {::},
issn = {1432-1416},
journal = {Journal of mathematical biology},
keywords = {Animals,Genetic Speciation,Genetics, Population,Humans,Likelihood Functions,Mathematical Concepts,Models, Genetic,Phylogeny,Population Density},
number = {1},
pages = {95--106},
pmid = {19283383},
title = {{Maximum tree: a consistent estimator of the species tree.}},
url = {http://www.ncbi.nlm.nih.gov/pubmed/19283383},
volume = {60},
year = {2010}
}

@article{Mossel2008,
abstract = {We introduce a simple computationally efficient algorithm for reconstructing phylogenies from multiple gene trees in the presence of incomplete lineage sorting, that is, when the topology of the gene trees may differ from that of the species tree. We show that our technique is statistically consistent under standard stochastic assumptions, that is, it returns the correct tree given sufficiently many unlinked loci. We also show that it can tolerate moderate estimation errors.},
author = {Mossel, Elchanan and Roch, Sebastien},
doi = {10.1109/TCBB.2008.66},
file = {:home/leo/articles/genome.evolutionary.theory/deepcoal/arxiv.0710.0262v2.pdf:pdf;:home/leo/articles/genome.evolutionary.theory/deepcoal/TransComputBiolBioinfo.2008.doi200866.pdf:pdf},
issn = {1557-9964},
journal = {IEEE/ACM transactions on computational biology and bioinformatics / IEEE, ACM},
number = {1},
pages = {166--71},
pmid = {20150678},
title = {{Incomplete lineage sorting: consistent phylogeny estimation from multiple loci.}},
url = {http://www.ncbi.nlm.nih.gov/pubmed/20150678},
volume = {7},
year = {2008}
}

@article{Liu2009,
abstract = {We review recent models to estimate phylogenetic trees under the multispecies coalescent. Although the distinction between gene trees and species trees has come to the fore of phylogenetics, only recently have methods been developed that explicitly estimate species trees. Of the several factors that can cause gene tree heterogeneity and discordance with the species tree, deep coalescence due to random genetic drift in branches of the species tree has been modeled most thoroughly. Bayesian approaches to estimating species trees utilizes two likelihood functions, one of which has been widely used in traditional phylogenetics and involves the model of nucleotide substitution, and the second of which is less familiar to phylogeneticists and involves the probability distribution of gene trees given a species tree. Other recent parametric and nonparametric methods for estimating species trees involve parsimony criteria, summary statistics, supertree and consensus methods. Species tree approaches are an appropriate goal for systematics, appear to work well in some cases where concatenation can be misleading, and suggest that sampling many independent loci will be paramount. Such methods can also be challenging to implement because of the complexity of the models and computational time. In addition, further elaboration of the simplest of coalescent models will be required to incorporate commonly known issues such as deviation from the molecular clock, gene flow and other genetic forces.},
author = {Liu, Liang and Yu, Lili and Kubatko, Laura Salter and Pearl, Dennis K and Edwards, Scott V},
doi = {10.1016/j.ympev.2009.05.033},
file = {:home/leo/articles/coalescent/MolPhylogEvol.2009.320doi05033.pdf:pdf},
issn = {1095-9513},
journal = {Molecular phylogenetics and evolution},
keywords = {Bayes Theorem,Evolution, Molecular,Genetic Speciation,Likelihood Functions,Models, Genetic,Phylogeny,Sequence Analysis, DNA,Sequence Analysis, DNA: methods,Statistics, Nonparametric},
number = {1},
pages = {320--8},
pmid = {19501178},
publisher = {Elsevier Inc.},
title = {{Coalescent methods for estimating phylogenetic trees.}},
url = {http://www.ncbi.nlm.nih.gov/pubmed/19501178},
volume = {53},
year = {2009}
}

@article{Kubatko2009b,
author = {Kubatko, Laura Salter and Carstens, Bryan C and Knowles, L Lacey},
doi = {10.1093/bioinformatics/btp079},
file = {:home/leo/articles/genome.evolutionary.theory/deepcoal/Bioinformatics.2009.971.pdf:pdf;:home/leo/articles/genome.evolutionary.theory/deepcoal/Bioinformatics.2009.971\_suppl.pdf:pdf},
issn = {1367-4811},
journal = {Bioinformatics},
keywords = {Algorithms,Computer Simulation,Genes,Likelihood Functions,Phylogeny,Software},
number = {7},
pages = {971--3},
pmid = {19211573},
title = {{STEM: species tree estimation using maximum likelihood for gene trees under coalescence.}},
url = {http://www.ncbi.nlm.nih.gov/pubmed/19211573},
volume = {25},
year = {2009}
}

% this is not on mendeley bibtex yet...
@article{Liu2010phybase,
author = {Liu, Liang and Yu, Lili},
title = {{Phybase: an R package for species tree analysis}},
journal = {Bioinformatics},
volume = {26},
number = {7},
pages = {962-963},
doi = {10.1093/bioinformatics/btq062},
year = {2010},
abstract = {Motivation: Phybase is an R package for phylogenetic analysis using species trees. It provides functions to read, write, manipulate, simulate, estimate, summarize and plot species trees, which contain not only the topology and branch lengths but also population sizes.  Availability: The Phybase package is available at the R repository. The manual and supporting materials including source code, sample R code and sample data files for the species tree analysis are available at http://stat.osu.edu/[~]liuliang/research/phybase.html  Contact: lliu@desu.edu
},
URL = {http://bioinformatics.oxfordjournals.org/cgi/content/abstract/26/7/962},
eprint = {http://bioinformatics.oxfordjournals.org/cgi/reprint/26/7/962.pdf}
}

@article{Liu2009d,
author = {Liu, L. and Edwards, Scott V},
doi = {10.1093/sysbio/syp034},
file = {:home/leo/articles/genome.evolutionary.theory/deepcoal/SystBiol.2009.452.pdf:pdf},
issn = {1063-5157},
journal = {Systematic Biology},
number = {4},
pages = {452--460},
title = {{Phylogenetic Analysis in the Anomaly Zone}},
url = {http://sysbio.oxfordjournals.org/cgi/doi/10.1093/sysbio/syp034},
volume = {58},
year = {2009}
}

@article{Liu2009e,
author = {Liu, L. and Yu, Lili and Pearl, Dennis K and Edwards, Scott V},
doi = {10.1093/sysbio/syp031},
file = {:home/leo/articles/genome.evolutionary.theory/deepcoal/SystBiol.2009.468syp031v1.pdf:pdf},
issn = {1063-5157},
journal = {Systematic Biology},
number = {5},
pages = {468--477},
title = {{Estimating Species Phylogenies Using Coalescence Times among Sequences}},
url = {http://sysbio.oxfordjournals.org/cgi/doi/10.1093/sysbio/syp031},
volume = {58},
year = {2009}
}

@article{McCormack2009b,
author = {McCormack, J. E. and Huang, H. and Knowles, L Lacey},
doi = {10.1093/sysbio/syp045},
file = {:home/leo/articles/genome.evolutionary.theory/deepcoal/SystBiol.2009.501syp045v1.pdf:pdf},
issn = {1063-5157},
journal = {Systematic Biology},
number = {5},
pages = {501--508},
title = {{Maximum Likelihood Estimates of Species Trees: How Accuracy of Phylogenetic Inference Depends upon the Divergence History and Sampling Design}},
url = {http://sysbio.oxfordjournals.org/cgi/doi/10.1093/sysbio/syp045},
volume = {58},
year = {2009}
}
% LCA mapping
@article{Zmasek2001,
author = {Zmasek, CM and Eddy, SR},
file = {:home/leo/articles/02-ImportedToMendeley/Bioinformatics.2001.821.pdf:pdf},
journal = {Bioinformatics},
number = {9},
pages = {821--828},
title = {{A simple algorithm to infer gene duplication and speciation events on a gene tree}},
url = {http://bioinformatics.oxfordjournals.org/content/17/9/821.short},
volume = {17},
year = {2001}
}
% cummulative number of losses below node (used on guenomu's source code)
@article{Ma2001,
author = {Ma, B and Li, Ming and Zhang, Louxin},
file = {:home/leo/articles/02-ImportedToMendeley/SIAM.2000.729.pdf:pdf},
journal = {SIAM Journal on Computing},
number = {3},
pages = {729--752},
title = {{From gene trees to species trees}},
url = {http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.95.7846\&amp;rep=rep1\&amp;type=pdf},
volume = {30},
year = {2001}
}
% cummulative number of losses over edges (described but not used on guenomu's source code)
@article{Tofigh2011,
abstract = {The incongruency between a gene tree and a corresponding species tree can be attributed to evolutionary events such as gene duplication and gene loss. This paper describes a combinatorial model where so-called DTL-scenarios are used to explain the differences between a gene tree and a corresponding species tree taking into account gene duplications, gene losses, and lateral gene transfers (also known as horizontal gene transfers). The reasonable biological constraint that a lateral gene transfer may only occur between contemporary species leads to the notion of acyclic DTL-scenarios. Parsimony methods are introduced by defining appropriate optimization problems. We show that finding most parsimonious acyclic DTL-scenarios is NP-hard. However, by dropping the condition of acyclicity, the problem becomes tractable, and we provide a dynamic programming algorithm as well as a fixed-parameter tractable algorithm for finding most parsimonious DTL-scenarios.},
author = {Tofigh, Ali and Hallett, Michael and Lagergren, Jens},
doi = {10.1109/TCBB.2010.14},
file = {:home/leo/articles/02-ImportedToMendeley/TofighAli-manuscript1-ACMTransComputBiolBioinfo.2010.14.pdf:pdf;:home/leo/articles/02-ImportedToMendeley/TofighAli-ACMTransComputBiolBioinfo.2011.517.pdf:pdf},
issn = {1557-9964},
journal = {IEEE/ACM Transactions on Computational Biology and Bioinformatics (TCBB)},
keywords = {Algorithms,DLT model,Evolution,Gene Duplication,Gene Transfer,Genetic,Horizontal,Horizontal: genetics,Models,Molecular,duploss,reconciliation},
mendeley-tags = {DLT model,duploss,reconciliation},
number = {2},
pages = {517--35},
pmid = {21233529},
publisher = {IEEE Computer Society Press},
title = {{Simultaneous identification of duplications and lateral gene transfers.}},
url = {http://dl.acm.org/citation.cfm?id=1936451 http://www.ncbi.nlm.nih.gov/pubmed/21233529},
volume = {8},
year = {2011}
}
% deepcoal calc previously used by guenomu. The comment was: "for "regular nodes" extra lineages above a node are equal to 
% |subtree rooted at node| - # coalescences below subtree (including root node of subtree) - 1 
@article{Than2009,
abstract = {In a 1997 seminal paper, W. Maddison proposed minimizing deep coalescences, or MDC, as an optimization criterion for inferring the species tree from a set of incongruent gene trees, assuming the incongruence is exclusively due to lineage sorting. In a subsequent paper, Maddison and Knowles provided and implemented a search heuristic for optimizing the MDC criterion, given a set of gene trees. However, the heuristic is not guaranteed to compute optimal solutions, and its hill-climbing search makes it slow in practice. In this paper, we provide two exact solutions to the problem of inferring the species tree from a set of gene trees under the MDC criterion. In other words, our solutions are guaranteed to find the tree that minimizes the total number of deep coalescences from a set of gene trees. One solution is based on a novel integer linear programming (ILP) formulation, and another is based on a simple dynamic programming (DP) approach. Powerful ILP solvers, such as CPLEX, make the first solution appealing, particularly for very large-scale instances of the problem, whereas the DP-based solution eliminates dependence on proprietary tools, and its simplicity makes it easy to integrate with other genomic events that may cause gene tree incongruence. Using the exact solutions, we analyze a data set of 106 loci from eight yeast species, a data set of 268 loci from eight Apicomplexan species, and several simulated data sets. We show that the MDC criterion provides very accurate estimates of the species tree topologies, and that our solutions are very fast, thus allowing for the accurate analysis of genome-scale data sets. Further, the efficiency of the solutions allow for quick exploration of sub-optimal solutions, which is important for a parsimony-based criterion such as MDC, as we show. We show that searching for the species tree in the compatibility graph of the clusters induced by the gene trees may be sufficient in practice, a finding that helps ameliorate the computational requirements of optimization solutions. Further, we study the statistical consistency and convergence rate of the MDC criterion, as well as its optimality in inferring the species tree. Finally, we show how our solutions can be used to identify potential horizontal gene transfer events that may have caused some of the incongruence in the data, thus augmenting Maddison's original framework. We have implemented our solutions in the PhyloNet software package, which is freely available at: http://bioinfo.cs.rice.edu/phylonet.},
author = {Than, Cuong and Nakhleh, Luay},
doi = {10.1371/journal.pcbi.1000501},
file = {:home/leo/articles/02-ImportedToMendeley/PLoScbi.2009.1000501.pdf:pdf;:home/leo/articles/02-ImportedToMendeley/PLoScbi.2009.1000501\_draft.pdf:pdf},
issn = {1553-7358},
journal = {PLoS computational biology},
keywords = {Animals,Apicomplexa,Apicomplexa: classification,Apicomplexa: genetics,Cluster Analysis,Computational Biology,Computational Biology: methods,Databases,Genetic,Genetic Speciation,Models,Phylogeny,Yeasts,Yeasts: classification,Yeasts: genetics,deep coalescence,phylogenomics},
mendeley-tags = {deep coalescence,phylogenomics},
month = sep,
number = {9},
pages = {e1000501},
pmid = {19749978},
title = {{Species tree inference by minimizing deep coalescences.}},
url = {http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2729383\&tool=pmcentrez\&rendertype=abstract http://www.ncbi.nlm.nih.gov/pubmed/19749978},
volume = {5},
year = {2009}
}
% obs for deepcoal calc previously used by guenomu (nowadays we use nloss - 2xndups + 2x|diff in size|). The notice was:
%" they mention the need for creating artificial nodes for duplicated species and removal of subtrees with species unrepresented in gene family. 
%  For deep coal, we look only at the species tree [1], where the mapping represents the coalescences. 
%  [1] the exception are maybe the species with more than one copy, but I didn't look at that. "
@article{Slowinski1997,
author = {Slowinski, J B and Knight, A and Rooney, A P},
file = {:home/leo/articles/02-ImportedToMendeley/MolPhylogEvol.1997.349.pdf:pdf},
journal = {Molecular Phylogenetics and Evolution},
number = {3},
pages = {349--362},
title = {{Inferring species trees from gene trees: a phylogenetic analysis of the Elapidae (Serpentes) based on the amino acid sequences of venom proteins}},
volume = {8},
year = {1997}
}
% effect of normalization constant
@article{Bryant2008,
abstract = {The Robinson-Foulds (RF) distance is by far the most widely used measure of dissimilarity between trees. Although the distribution of these distances has been investigated for twenty years, an algorithm that is explicitly polynomial time has yet to be described for computing this distribution (which is also the distribution of trees around a given tree under the popular Robinson-Foulds metric). In this paper we derive a polynomial-time algorithm for this distribution. We show how the distribution can be approximated by a Poisson distribution determined by the proportion of leaves that lie in `cherries' of the given tree. We also describe how our results can be used to derive normalization constants that are required in a recently-proposed maximum likelihood approach to supertree construction.},
archivePrefix = {arXiv},
arxivId = {0810.0868},
author = {Bryant, David and Steel, Mike},
eprint = {0810.0868},
file = {:home/leo/articles/02-ImportedToMendeley/arxiv.0810.0868-ACMTransComputBiolBioinfo.2009.420.pdf:pdf;:home/leo/articles/02-ImportedToMendeley/ACMTransComputBiolBioinfo.2009.420.pdf:pdf},
journal = {Discrete Mathematics},
keywords = {Populations and Evolution,Quantitative Methods},
month = oct,
pages = {420--426},
volume = {6},
number = {3},
doi = {10.1109/TCBB.2009.32},
title = {{Computing the Distribution of a Tree Metric}},
url = {http://arxiv.org/abs/0810.0868},
year = {2008}
}

@Article{Chaudhary2013,
  Title                    = {{Inferring species trees from incongruent multi-copy gene trees using the Robinson-Foulds distance}},
  Author                   = {Chaudhary, Ruchi and Burleigh, J. Gordon and Fern\'andez-Baca, David},
  Journal                  = {Algorithms for Molecular Biology},
  Year                     = {2013},
  Number                   = {1},
  Pages                    = {28},
  Volume                   = {8},

  Doi                      = {10.1186/1748-7188-8-28},
  File                     = {:home/leo/articles/021/AlgoMolBiol.2013.8-28.pdf:PDF},
  Keywords                 = {00.read,SPR distance, 00.distance},
  Owner                    = {leo},
  Publisher                = {Springer (Biomed Central Ltd.)},
  Timestamp                = {20131103.1250},
  Url                      = {http://dx.doi.org/10.1186/1748-7188-8-28}
}

@Article{Helmkamp2012,
  Title                    = {Improvements to a class of distance matrix methods for inferring species trees from gene trees.},
  Author                   = {Helmkamp, Laura J and Jewett, Ethan M and Rosenberg, Noah A},
  Journal                  = {J. Comput. Biol.},
  Year                     = {2012},
  Number                   = {6},
  Pages                    = {632-49},
  Volume                   = {19},
  Doi                      = {10.1089/cmb.2012.0042},
  File                     = {file:///home/leo/articles/02-ImportedToMendeley/JComputBiol.2012.632.pdf:pdf},
  Owner                    = {leo},
  Pmid                     = {22697239},
  Timestamp                = {2013.08.18}
}

@Article{Chaudhary2010,
  Title                    = {{iGTP: a software package for large-scale gene tree parsimony analysis.}},
  Author                   = {Chaudhary, Ruchi and Bansal, Mukul S and Wehe, Andr\'{e} and Fern\'{a}ndez-Baca, David and Eulenstein, Oliver},
  Journal                  = {BMC bioinformatics},
  Year                     = {2010},
  Month                    = jan,
  Number                   = {1},
  Pages                    = {574},
  Volume                   = {11},
  Abstract                 = {The ever-increasing wealth of genomic sequence information provides an unprecedented opportunity for large-scale phylogenetic analysis. However, species phylogeny inference is obfuscated by incongruence among gene trees due to evolutionary events such as gene duplication and loss, incomplete lineage sorting (deep coalescence), and horizontal gene transfer. Gene tree parsimony (GTP) addresses this issue by seeking a species tree that requires the minimum number of evolutionary events to reconcile a given set of incongruent gene trees. Despite its promise, the use of gene tree parsimony has been limited by the fact that existing software is either not fast enough to tackle large data sets or is restricted in the range of evolutionary events it can handle.},
  Doi                      = {10.1186/1471-2105-11-574},
  File                     = {:home/leo/articles/02-ImportedToMendeley/BMCBioinfo.2010.11-574.pdf:pdf},
  ISSN                     = {1471-2105},
  Keywords                 = {Algorithms,Databases, Genetic,Evolution, Molecular,Gene Duplication,Genome,Genomics,Genomics: methods,Phylogeny,Software},
  Pmid                     = {21092314},
  Publisher                = {BioMed Central Ltd},
  Url                      = {http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3002902\&tool=pmcentrez\&rendertype=abstract}
}

@Article{Mallo2013,
  Title   = {SimPhy: Comprehensive simulation of gene, locus and species trees at the genome-wide level},
  Author  = {Mallo, Diego and de Oliveira Martins, Leonardo and Posada, David},
  Journal = {unpublished},
  Year    = {2014},
  Url     = {https://bitbucket.org/dmallo/stidsim}
}

@article{Sukumaran2010,
 title     ={DendroPy: a Python library for phylogenetic computing},
 author    ={Sukumaran, Jeet and Holder, Mark T},
 journal   ={Bioinformatics},
 volume    ={26},
 number    ={12},
 pages     ={1569--1571},
 year      ={2010},
 Doi       ={10.1093/bioinformatics/btq228},
 Url       ={http://bioinformatics.oxfordjournals.org/content/26/12/1569.full},
 publisher ={Oxford Univ Press}
}

@Article{Leache2011,
  Title                    = {{The accuracy of species tree estimation under simulation: a comparison of methods.}},
  Author                   = {Leach\'{e}, Adam D and Rannala, Bruce},
  Journal                  = {Systematic biology},
  Year                     = {2011},
  Month                    = mar,
  Number                   = {2},
  Pages                    = {126--37},
  Volume                   = {60},
  Doi                      = {10.1093/sysbio/syq073},
  File                     = {:home/leo/articles/02-ImportedToMendeley/SystBiol.2011.126.pdf:pdf},
  ISSN                     = {1076-836X},
  Keywords                 = {Bayes Theorem,Classification,Classification: methods,Computer Simulation,Genetic Speciation,Phylogeny,Plants,Plants: classification,Plants: genetics},
  Pmid                     = {21088009},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/21088009}
}

@article{Ane2007,
  title={Bayesian estimation of concordance among gene trees},
  author={An{\'e}, C{\'e}cile and Larget, Bret and Baum, David A and Smith, Stacey D and Rokas, Antonis},
  journal={Molecular Biology and Evolution},
  volume={24},
  number={2},
  pages={412--426},
  year={2007},
  publisher={SMBE}
}

@Article{DeGiorgio2013,
  Title                    = {Robustness to Divergence Time Underestimation When Inferring Species Trees from Estimated Gene Trees},
  Author                   = {DeGiorgio, M. and Degnan, J. H.},
  Journal                  = {Systematic Biology},
  Year                     = {2013},

  Month                    = {Aug},

  Doi                      = {10.1093/sysbio/syt059},
  File                     = {:home/leo/articles/021/SystBiol.2013.syt059.pdf:PDF},
  Owner                    = {leo},
  Publisher                = {Oxford University Press},
  Timestamp                = {2013.10.09-14:30},
  Url                      = {http://dx.doi.org/10.1093/sysbio/syt059}
}

@Article{Steel2008,
  Title                    = {{Maximum likelihood supertrees.}},
  Author                   = {Steel, Mike and Rodrigo, Allen},
  Journal                  = {Systematic biology},
  Year                     = {2008},

  Month                    = apr,
  Number                   = {2},
  Pages                    = {243--50},
  Volume                   = {57},

  Abstract                 = {We analyze a maximum likelihood approach for combining phylogenetic trees into a larger "supertree." This is based on a simple exponential model of phylogenetic error, which ensures that ML supertrees have a simple combinatorial description (as a median tree, minimizing a weighted sum of distances to the input trees). We show that this approach to ML supertree reconstruction is statistically consistent (it converges on the true species supertree as more input trees are combined), in contrast to the widely used MRP method, which we show can be statistically inconsistent under the exponential error model. We also show that this statistical consistency extends to an ML approach for constructing species supertrees from gene trees. In this setting, incomplete lineage sorting (due to coalescence rates of homologous genes being lower than speciation rates) has been shown to lead to gene trees that are frequently different from species trees, and this can confound efforts to reconstruct the species phylogeny correctly.},
  Archiveprefix            = {arXiv},
  Arxivid                  = {0708.2124},
  Doi                      = {10.1080/10635150802033014},
  Eprint                   = {0708.2124},
  File                     = {:home/leo/articles/02-ImportedToMendeley/SystBiol.2008.243.pdf:pdf},
  ISSN                     = {1063-5157},
  Keywords                 = {Computer Simulation,Likelihood Functions,Models, Genetic,Phylogeny},
  Pmid                     = {18398769},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/18398769}
}

@Article{Cotton2009,
  Title                    = {{Supertrees join the mainstream of phylogenetics}},
  Author                   = {Cotton, James A and Wilkinson, M.},
  Journal                  = {Trends in ecology \& evolution},
  Year                     = {2009},
  Number                   = {1},
  Pages                    = {1--3},
  Volume                   = {24},

  File                     = {:home/leo/articles/02-ImportedToMendeley/TrendsEcolEvol.2008.doi08006-2009.1.pdf:pdf},
  Keywords                 = {bayesian,genome evolution,genomics,phylogenomics},
  Mendeley-tags            = {genome evolution,genomics,phylogenomics},
  Publisher                = {Elsevier},
  Url                      = {http://linkinghub.elsevier.com/retrieve/pii/S0169534708003261}
}

@Article{Akerborg2009,
  Title                    = {{Simultaneous Bayesian gene tree reconstruction and reconciliation analysis.}},
  Author                   = {Akerborg, Orjan and Sennblad, Bengt and Arvestad, Lars and Lagergren, Jens},
  Journal                  = {Proceedings of the National Academy of Sciences of the United States of America},
  Year                     = {2009},

  Month                    = apr,
  Number                   = {14},
  Pages                    = {5714--9},
  Volume                   = {106},

  Abstract                 = {We present GSR, a probabilistic model integrating gene duplication, sequence evolution, and a relaxed molecular clock for substitution rates, that enables genomewide analysis of gene families. The gene duplication and loss process is a major cause for incongruence between gene and species tree, and deterministic methods have been developed to explain such differences through tree reconciliations. Although probabilistic methods for phylogenetic inference have been around for decades, probabilistic reconciliation methods are far less established. Based on our model, we have implemented a Bayesian analysis tool, PrIME-GSR, for gene tree inference that takes a known species tree into account. Our implementation is sound and we demonstrate its utility for genomewide gene-family analysis by applying it to recently presented yeast data. We validate PrIME-GSR by comparing with previous analyses of these data that take advantage of gene order information. In a case study we apply our method to the ADH gene family and are able to draw biologically relevant conclusions concerning gene duplications creating key yeast phenotypes. On a higher level this shows the biological relevance of our method. The obtained results demonstrate the value of a relaxed molecular clock. Our good performance will extend to species where gene order conservation is insufficient.},
  Doi                      = {10.1073/pnas.0806251106},
  File                     = {:home/leo/articles/02-ImportedToMendeley/PNAS.2009.5714.pdf:pdf;:home/leo/articles/02-ImportedToMendeley/PNAS.2009.5714\_Suppl.pdf:pdf},
  ISSN                     = {1091-6490},
  Keywords                 = {Bayes Theorem,Fungal,Gene Duplication,Genetic,Genome,Kinetics,Models,Mutation,Phylogeny,Yeasts,Yeasts: genetics,duploss,probabilistic,reconciliation},
  Mendeley-tags            = {duploss,probabilistic,reconciliation},
  Pmid                     = {19299507},
  Url                      = {http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2667006\&tool=pmcentrez\&rendertype=abstract}
}

@Article{Boussau2013,
  Title                    = {Genome-scale coestimation of species and gene trees.},
  Author                   = {Boussau, Bastien and Sz{\"{o}}llosi, Gergely J and Duret, Laurent and Gouy, Manolo and Tannier, Eric and Daubin, Vincent},
  Journal                  = {Genome Res.},
  Year                     = {2013},
  Number                   = {2},
  Pages                    = {323-30},
  Volume                   = {23},

  Abstract                 = {Comparisons of gene trees and species trees are key to understanding major processes of genome evolution such as gene duplication and loss. Because current methods to reconstruct phylogenies fail to model the two-way dependency between gene trees and the species tree, they often misrepresent gene and species histories. We present a new probabilistic model to jointly infer rooted species and gene trees for dozens of genomes and thousands of gene families. We use simulations to show that this method accurately infers the species tree and gene trees, is robust to misspecification of the models of sequence and gene family evolution, and provides a precise historic record of gene duplications and losses throughout genome evolution. We simultaneously reconstruct the history of mammalian species and their genes based on 36 completely sequenced genomes, and use the reconstructed gene trees to infer the gene content and organization of ancestral mammalian genomes. We show that our method yields a more accurate picture of ancestral genomes than the trees available  in the authoritative database Ensembl.},
  Doi                      = {10.1101/gr.141978.112},
  File                     = {file:///home/leo/articles/02-ImportedToMendeley/GenomeRes.2013.323.pdf:pdf},
  Owner                    = {leo},
  Pmid                     = {23132911},
  Timestamp                = {2013.08.18}
}

@Article{Heled2010,
  Title                    = {{Bayesian inference of species trees from multilocus data.}},
  Author                   = {Heled, Joseph and Drummond, Alexei J},
  Journal                  = {Molecular biology and evolution},
  Year                     = {2010},

  Month                    = mar,
  Number                   = {3},
  Pages                    = {570--80},
  Volume                   = {27},

  Abstract                 = {Until recently, it has been common practice for a phylogenetic analysis to use a single gene sequence from a single individual organism as a proxy for an entire species. With technological advances, it is now becoming more common to collect data sets containing multiple gene loci and multiple individuals per species. These data sets often reveal the need to directly model intraspecies polymorphism and incomplete lineage sorting in phylogenetic estimation procedures. For a single species, coalescent theory is widely used in contemporary population genetics to model intraspecific gene trees. Here, we present a Bayesian Markov chain Monte Carlo method for the multispecies coalescent. Our method coestimates multiple gene trees embedded in a shared species tree along with the effective population size of both extant and ancestral species. The inference is made possible by multilocus data from multiple individuals per species. Using a multiindividual data set and a series of simulations of rapid species radiations, we demonstrate the efficacy of our new method. These simulations give some insight into the behavior of the method as a function of sampled individuals, sampled loci, and sequence length. Finally, we compare our new method to both an existing method (BEST 2.2) with similar goals and the supermatrix (concatenation) method. We demonstrate that both BEST and our method have much better estimation accuracy for species tree topology than concatenation, and our method outperforms BEST in divergence time and population size estimation.},
  Doi                      = {10.1093/molbev/msp274},
  File                     = {:home/leo/articles/02-ImportedToMendeley/MolBiolEvol.2010.570.pdf:PDF;:home/leo/articles/02-ImportedToMendeley/MolBiolEvol.2010.570-draft_msp274v2.pdf:PDF},
  ISSN                     = {1537-1719},
  Keywords                 = {Algorithms,Bayes Theorem,Computer Simulation,Evolution,Genes,Genetic,Models,Molecular,Monte Carlo Method,Phylogeny,Species Specificity,bayesian inference,censored,deep coalescence,gene trees,molecular systematics,multispecies coalescent,phylogenomics,species trees},
  Mendeley-tags            = {deep coalescence,phylogenomics},
  Pmid                     = {19906793},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/19906793 http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2822290\&tool=pmcentrez\&rendertype=abstract}
}

@Article{Edwards2007,
  Title                    = {{High-resolution species trees without concatenation.}},
  Author                   = {Edwards, Scott V and Liu, Liang and Pearl, Dennis K},
  Journal                  = {Proceedings of the National Academy of Sciences of the United States of America},
  Year                     = {2007},

  Month                    = apr,
  Number                   = {14},
  Pages                    = {5936--41},
  Volume                   = {104},

  Abstract                 = {The vast majority of phylogenetic models focus on resolution of gene trees, despite the fact that phylogenies of species in which gene trees are embedded are of primary interest. We analyze a Bayesian model for estimating species trees that accounts for the stochastic variation expected for gene trees from multiple unlinked loci sampled from a single species history after a coalescent process. Application of the model to a 106-gene data set from yeast shows that the set of gene trees recovered by statistically acknowledging the shared but unknown species tree from which gene trees are sampled is much reduced compared with treating the history of each locus independently of an overarching species tree. The analysis also yields a concentrated posterior distribution of the yeast species tree whose mode is congruent with the concatenated gene tree but can do so with less than half the loci required by the concatenation method. Using simulations, we show that, with large numbers of loci, highly resolved species trees can be estimated under conditions in which concatenation of sequence data will positively mislead phylogeny, and when the proportion of gene trees matching the species tree is <10\%. However, when gene tree/species tree congruence is high, species trees can be resolved with just two or three loci. These results make accessible an alternative paradigm for combining data in phylogenomics that focuses attention on the singularity of species histories and away from the idiosyncrasies and multiplicities of individual gene histories.},
  Doi                      = {10.1073/pnas.0607004104},
  File                     = {:home/leo/articles/02-ImportedToMendeley/PNAS.2007.5936.pdf:pdf;:home/leo/articles/02-ImportedToMendeley/PNAS.2007.5936\_suppl.pdf:pdf},
  ISSN                     = {0027-8424},
  Keywords                 = {Alleles,Base Composition,Base Sequence,Bayes Theorem,Concatenated,DNA,Fungal,Fungal: genetics,Genetic,Genomics,Likelihood Functions,Markov Chains,Models,Molecular Sequence Data,Monte Carlo Method,Phylogeny,Saccharomyces cerevisiae,Saccharomyces cerevisiae: genetics,Sequence Analysis,Species Specificity,Statistical,Stochastic Processes,deep coalescence,phylogenomics},
  Mendeley-tags            = {deep coalescence,phylogenomics},
  Pmid                     = {17392434},
  Url                      = {http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=1851595\&tool=pmcentrez\&rendertype=abstract}
}

@Article{Schreiber2013,
  Title                    = {TreeFam v9: a new website, more species and orthology-on-the-fly},
  Author                   = {Schreiber, F. and Patricio, M. and Muffato, M. and Pignatelli, M. and Bateman, A.},
  Journal                  = {Nucleic Acids Research},
  Year                     = {2014},
  Volume = {42},
  Number = {D1},
  Pages = {D922--D925},
  Month                    = {Nov},

  Doi                      = {10.1093/nar/gkt1055},
  File                     = {:home/leo/articles/021/NuclAcidsRes.2013.gkt1055.pdf:PDF},
  Owner                    = {leo},
  Publisher                = {Oxford University Press},
  Timestamp                = {20131106.1809},
  Url                      = {http://dx.doi.org/10.1093/nar/gkt1055}
}

@InCollection{Ronquist2004,
  Title                    = {Bayesian supertrees},
  Author                   = {Ronquist, Fredrik and Huelsenbeck, John P. and Britton, Tom},
  Booktitle                = {Phylogenetic supertrees},
  Publisher                = {Springer},
  Year                     = {2004},
  Editor                   = {Olaf R. P. Bininda-Emonds},
  Pages                    = {193-224},
  File                     = {file:///home/leo/articles/02-ImportedToMendeley/BayesianSupertrees.2004.193.pdf:pdf},
  Owner                    = {leo},
  Timestamp                = {2013.07.01}
}

@Article{Sjostrand2012,
  Title                    = {{DLRS: gene tree evolution in light of a species tree}},
  Author                   = {Sjostrand, J. and Sennblad, B. and Arvestad, L. and Lagergren, J.},
  Journal                  = {Bioinformatics},
  Year                     = {2012},

  Month                    = {Nov},
  Number                   = {22},
  Pages                    = {2994-2995},
  Volume                   = {28},

  Doi                      = {10.1093/bioinformatics/bts548},
  File                     = {:home/leo/articles/021/Bioinformatics.2012.2994.pdf:PDF},
  Keywords                 = {00.reconciliation},
  Owner                    = {leo},
  Publisher                = {Oxford University Press},
  Timestamp                = {20131120.1413},
  Url                      = {http://dx.doi.org/10.1093/bioinformatics/bts548}
}


@Article{Degiorgio2010,
  Title                    = {{Fast and consistent estimation of species trees using supermatrix rooted triples.}},
  Author                   = {Degiorgio, Michael and Degnan, James H},
  Journal                  = {Molecular biology and evolution},
  Year                     = {2010},
  Number                   = {3},
  Pages                    = {552--69},
  Volume                   = {27},

  Abstract                 = {Concatenated sequence alignments are often used to infer species-level relationships. Previous studies have shown that analysis of concatenated data using maximum likelihood (ML) can produce misleading results when loci have differing gene tree topologies due to incomplete lineage sorting. Here, we develop a polynomial time method that utilizes the modified mincut supertree algorithm to construct an estimated species tree from inferred rooted triples of concatenated alignments. We term this method SuperMatrix Rooted Triple (SMRT) and use the notation SMRT-ML when rooted triples are inferred by ML. We use simulations to investigate the performance of SMRT-ML under Jukes-Cantor and general time-reversible substitution models for four- and five-taxon species trees and also apply the method to an empirical data set of yeast genes. We find that SMRT-ML converges to the correct species tree in many cases in which ML on the full concatenated data set fails to do so. SMRT-ML can be conservative in that its output tree is often partially unresolved for problematic clades. We show analytically that when the species tree is clocklike and mutations occur under the Cavender-Farris-Neyman substitution model, as the number of genes increases, SMRT-ML is increasingly likely to infer the correct species tree even when the most likely gene tree does not match the species tree. SMRT-ML is therefore a computationally efficient and statistically consistent estimator of the species tree when gene trees are distributed according to the multispecies coalescent model.},
  File                     = {:home/leo/articles/02-ImportedToMendeley/MolBiolEvol.2010.552.pdf:pdf},
  Keywords                 = {anomalous gene tree,anomaly zone,lineage sorting,phylogenetics,phylogenomics,statistical consistency},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/19833741}
}

@Article{Whidden2013,
  Title                    = {Supertrees based on the subtree prune-and-regraft distance},
  Author                   = {Whidden, Chris and Zeh, Norbert and Beiko, Robert G.},
  Journal                  = {PeerJ PrePrints},
  Year                     = {2013},

  Month                    = {5},
  Pages                    = {18},
  Volume                   = {1},

  Abstract                 = {Supertree methods reconcile a set of phylogenetic trees into a single structure that is often interpreted as a branching history of species. A key challenge is combining conflicting evolutionary histories that are due to artifacts of phylogenetic reconstruction and phenomena such as lateral gene transfer (LGT). Although they often work well in practice, existing supertree approaches use optimality criteria that do not reflect underlying processes, have known biases and may be unduly influenced by LGT. We present the first method to construct supertrees by using the subtree prune-and-regraft (SPR) distance as an optimality criterion. Although calculating the rooted SPR distance between a pair of trees is NP-hard, our new maximum agreement forest-based methods can reconcile trees with hundreds of taxa and > 50 transfers in fractions of a second, which enables repeated calculations during the course of an iterative search. Our approach can accommodate trees in which uncertain relationships have been collapsed to multifurcating nodes. Using a series of simulated benchmark datasets, we show that SPR supertrees are more similar to correct species histories under plausible rates of LGT than supertrees based on parsimony or Robinson-Foulds distance criteria. We successfully constructed an SPR supertree from a phylogenomic dataset of 40,631 gene trees that covered 244 genomes representing several major bacterial phyla. Our SPR-based approach also allowed direct inference of highways of gene transfer between bacterial classes and genera; a small number of these highways connect genera in different phyla and can highlight specific genes implicated in long-distance LGT.},
  Doi                      = {10.7287/peerj.preprints.18v1},
  File                     = {file:///home/leo/articles/02-ImportedToMendeley/peerj-preprints-18.2013.pdf:pdf},
  ISSN                     = {2167-9843},
  Keywords                 = {prokaryotic phylogeny, supertrees, matrix representation with parsimony, subtree prune-and-regraft, phylogenomics, lateral gene transfer, Robinson-Foulds},
  Owner                    = {leo},
  Timestamp                = {2013.06.03_12:51},
  Url                      = {http://dx.doi.org/10.7287/peerj.preprints.18v1}
}

@Article{Rannala2003,
  Title                    = {Bayes estimation of species divergence times and ancestral population sizes using DNA sequences from multiple loci.},
  Author                   = {Bruce Rannala and Ziheng Yang},
  Journal                  = {Genetics},
  Year                     = {2003},

  Month                    = {Aug},
  Number                   = {4},
  Pages                    = {1645--1656},
  Volume                   = {164},

  Abstract                 = {The effective population sizes of ancestral as well as modern species are important parameters in models of population genetics and human evolution. The commonly used method for estimating ancestral population sizes, based on counting mismatches between the species tree and the inferred gene trees, is highly biased as it ignores uncertainties in gene tree reconstruction. In this article, we develop a Bayes method for simultaneous estimation of the species divergence times and current and ancestral population sizes. The method uses DNA sequence data from multiple loci and extracts information about conflicts among gene tree topologies and coalescent times to estimate ancestral population sizes. The topology of the species tree is assumed known. A Markov chain Monte Carlo algorithm is implemented to integrate over uncertain gene trees and branch lengths (or coalescence times) at each locus as well as species divergence times. The method can handle any species tree and allows different numbers of sequences at different loci. We apply the method to published noncoding DNA sequences from the human and the great apes. There are strong correlations between posterior estimates of speciation times and ancestral population sizes. With the use of an informative prior for the human-chimpanzee divergence date, the population size of the common ancestor of the two species is estimated to be approximately 20,000, with a 95\% credibility interval (8000, 40,000). Our estimates, however, are affected by model assumptions as well as data quality. We suggest that reliable estimates have yet to await more data and more realistic models.},
  File                     = {:home/leo/articles/02-ImportedToMendeley/Genetics.2003.1645.pdf:pdf},
  Institution              = {Department of Medical Genetics, University of Alberta, Edmonton, Alberta T6G 2H7, Canada.},
  Keywords                 = {Algorithms; Alleles; Animals; Base Sequence; Bayes Theorem; Evolution, Molecular; Genetics, Population; Hominidae, genetics; Humans; Markov Chains; Models, Statistical; Monte Carlo Method; Pan troglodytes, genetics; Population Density; Reproducibility of Results; Species Specificity; Time Factors},
  Medline-pst              = {ppublish},
  Owner                    = {leo},
  Pmid                     = {12930768},
  Timestamp                = {2013.05.03}
}

@Article{Ronquist2012,
  Title                    = {{MrBayes 3.2: Efficient Bayesian Phylogenetic Inference and Model Choice Across a Large Model Space.}},
  Author                   = {Ronquist, Fredrik and Teslenko, Maxim and van der Mark, Paul and Ayres, Daniel L and
          Darling, Aaron and H\"{o}hna, Sebastian and Larget, Bret and Liu, Liang and Suchard, Marc A and Huelsenbeck, John P},
  Journal                  = {Systematic biology},
  Year                     = {2012},

  Month                    = feb,
  Pages                    = {1--14},

  Abstract                 = {Since its introduction in 2001, MrBayes has grown in popularity as a software package for Bayesian phylogenetic inference using MCMC methods. With this note we announce the release of version 3.2, a major upgrade to the latest official release presented in 2003. The new version provides convergence diagnostics and allows multiple analyses to be run in parallel with convergence progress monitored on the fly. The introduction of new proposals and automatic optimization of tuning parameters has improved convergence for many problems. The new version also sports significantly faster likelihood calculations through streaming single-instruction-multiple-data extensions (SSE) and support of the BEAGLE library, allowing likelihood calculations to be delegated to graphics processing units (GPUs) on compatible hardware. Speedup factors range from around 2 with SSE code to more than 50 with BEAGLE for codon problems. Check-pointing across all models allows long runs to be completed even when an analysis is prematurely terminated. New models include relaxed clocks, dating, model averaging across time-reversible substitution models, and support for hard, negative and partial (backbone) tree constraints. Inference of species trees from gene trees is supported by full incorporation of the Bayesian estimation of species trees (BEST) algorithms. Marginal model likelihoods for Bayes factor tests can be estimated accurately across the entire model space using the stepping-stone method. The new version provides more output options than previously, including samples of ancestral states, site rates, site d(N)/d(S) rations, branch rates, and node dates. A wide range of statistics on tree parameters can also be output for visualization in FigTree and compatible software.},
  Doi                      = {10.1093/sysbio/sys029},
  File                     = {:home/leo/articles/02-ImportedToMendeley/SystBiol.2012.sys029.pdf:pdf},
  ISSN                     = {1076-836X},
  Pmid                     = {22357727},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/22357727}
}

@Article{Arvestad2003,
  Title                    = {{Bayesian gene/species tree reconciliation and orthology analysis using MCMC}},
  Author                   = {Arvestad, L. and a.-C. Berglund and Lagergren, J. and Sennblad, B.},
  Journal                  = {Bioinformatics},
  Year                     = {2003},

  Month                    = jul,
  Number                   = {Suppl 1},
  Pages                    = {i7--i15},
  Volume                   = {19},

  Doi                      = {10.1093/bioinformatics/btg1000},
  File                     = {:home/leo/articles/02-ImportedToMendeley/Bioinformatics.2003.i7.pdf:pdf},
  ISSN                     = {1367-4803},
  Keywords                 = {duploss,probabilistic,reconciliation},
  Mendeley-tags            = {duploss,probabilistic,reconciliation},
  Url                      = {http://bioinformatics.oxfordjournals.org/cgi/doi/10.1093/bioinformatics/btg1000}
}

@Article{Rasmussen2011,
  Title                    = {{A Bayesian approach for fast and accurate gene tree reconstruction.}},
  Author                   = {Rasmussen, Matthew D and Kellis, Manolis},
  Journal                  = {Molecular biology and evolution},
  Year                     = {2011},

  Month                    = jan,
  Number                   = {1},
  Pages                    = {273--90},
  Volume                   = {28},

  Abstract                 = {Recent sequencing and computing advances have enabled phylogenetic analyses to expand to both entire genomes and large clades, thus requiring more efficient and accurate methods designed specifically for the phylogenomic context. Here, we present SPIMAP, an efficient Bayesian method for reconstructing gene trees in the presence of a known species tree. We observe many improvements in reconstruction accuracy, achieved by modeling multiple aspects of evolution, including gene duplication and loss (DL) rates, speciation times, and correlated substitution rate variation across both species and loci. We have implemented and applied this method on two clades of fully sequenced species, 12 Drosophila and 16 fungal genomes as well as simulated phylogenies and find dramatic improvements in reconstruction accuracy as compared with the most popular existing methods, including those that take the species tree into account. We find that reconstruction inaccuracies of traditional phylogenetic methods overestimate the number of DL events by as much as 2-3-fold, whereas our method achieves significantly higher accuracy. We feel that the results and methods presented here will have many important implications for future investigations of gene evolution.},
  Doi                      = {10.1093/molbev/msq189},
  File                     = {:home/leo/articles/02-ImportedToMendeley/MolBiolEvol.2011.273.pdf:pdf},
  ISSN                     = {1537-1719},
  Keywords                 = {Algorithms,Animals,Bayes Theorem,Computational Biology,Computational Biology: methods,Evolution, Molecular,Gene Duplication,Genes, Fungal,Genes, Insect,Humans,Models, Genetic,Phylogeny},
  Pmid                     = {20660489},
  Url                      = {http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3002250\&tool=pmcentrez\&rendertype=abstract}
}

@Article{Szollosi2013,
  Title                    = {Efficient Exploration of the Space of Reconciled Gene Trees},
  Author                   = {Szollosi, G. J. and Rosikiewicz, W. and Boussau, B. and Tannier, E. and Daubin, V.},
  Journal                  = {Systematic Biology},
  Year                     = {2013},

  Month                    = {Oct},
  Number                   = {6},
  Pages                    = {901-912},
  Volume                   = {62},

  Doi                      = {10.1093/sysbio/syt054},
  File                     = {:home/leo/articles/021/SystBiol.2013.901.pdf:PDF;:home/leo/articles/02-ImportedToMendeley/arxiv.1306.2167v1.pdf:PDF},
  Publisher                = {Oxford University Press},
  Url                      = {http://dx.doi.org/10.1093/sysbio/syt054}
}

@Article{Rasmussen2012,
  Title                    = {{Unified modeling of gene duplication, loss, and coalescence using a locus tree}},
  Author                   = {Rasmussen, Matthew D and Kellis, Manolis},
  Journal                  = {Genome Research},
  Year                     = {2012},

  Month                    = jan,
  Number                   = {4},
  Pages                    = {755--65},
  Volume                   = {22},

  Abstract                 = {Gene phylogenies provide a rich source of information about the way evolution shapes genomes, populations, and phenotypes. In addition to substitutions, evolutionary events such as gene duplication and loss (as well as horizontal transfer) play a major role in gene evolution, and many phylogenetic models have been developed in order to reconstruct and study these events. However, these models typically make the simplifying assumption that population-related effects such as incomplete lineage sorting (ILS) are negligible. While this assumption may have been reasonable in some settings, it has become increasingly problematic as increased genome sequencing has led to denser phylogenies, where effects such as ILS are more prominent. To address this challenge, we present a new probabilistic model, DLCoal, that defines gene duplication and loss in a population setting, such that coalescence and ILS can be directly addressed. Interestingly, this model implies that in addition to the usual gene tree and species tree, there exists a third tree, the locus tree, which will likely have many applications. Using this model, we develop the first general reconciliation method that accurately infers gene duplications and losses in the presence of ILS, and we show its improved inference of orthologs, paralogs, duplications, and losses for a variety of clades, including flies, fungi, and primates. Also, our simulations show that gene duplications increase the frequency of ILS, further illustrating the importance of a joint model. Going forward, we believe that this unified model can offer insights to questions in both phylogenetics and population genetics.},
  Doi                      = {10.1101/gr.123901.111},
  File                     = {:home/leo/articles/02-ImportedToMendeley/GenomeRes.2012.doi111.pdf:pdf;:home/leo/articles/02-ImportedToMendeley/GenomeRes.2012.doi111\_suppl.pdf:pdf},
  ISSN                     = {1088-9051},
  Keywords                 = {00.read,Algorithms,Animals,Evolution,Gene Deletion,Gene Duplication,Gene Transfer,Genetic,Genetic Loci,Genetic Loci: genetics,Genome,Genome: genetics,Horizontal,Humans,Models,Molecular,Mutation,Phylogeny,Species Specificity,Statistical,Yeasts,Yeasts: classification,Yeasts: genetics,coalescent,deep coalescence,duploss,gene family,gene trees,phylogenomics,reconciliation,species tree},
  Mendeley-tags            = {00.read,coalescent,deep coalescence,duploss,gene family,gene trees,phylogenomics,reconciliation,species tree},
  Owner                    = {leo},
  Pmid                     = {22271778},
  Timestamp                = {2013.10.14-1518},
  Url                      = {http://genome.cshlp.org/cgi/doi/10.1101/gr.123901.111 http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3317157\&tool=pmcentrez\&rendertype=abstract}
}

@Article{Meng2009,
  Title                    = {{Detecting hybrid speciation in the presence of incomplete lineage sorting using gene tree incongruence: a model.}},
  Author                   = {Meng, Chen and Kubatko, Laura Salter},
  Journal                  = {Theoretical population biology},
  Year                     = {2009},
  Number                   = {1},
  Pages                    = {35--45},
  Volume                   = {75},

  Abstract                 = {The application of phylogenetic inference methods, to data for a set of independent genes sampled randomly throughout the genome, often results in substantial incongruence in the single-gene phylogenetic estimates. Among the processes known to produce discord between single-gene phylogenies, two of the best studied in a phylogenetic context are hybridization and incomplete lineage sorting. Much recent attention has focused on the development of methods for estimating species phylogenies in the presence of incomplete lineage sorting, but phylogenetic models that allow for hybridization have been more limited. Here we propose a model that allows incongruence in single-gene phylogenies to be due to both hybridization and incomplete lineage sorting, with the goal of determining the contribution of hybridization to observed gene tree incongruence in the presence of incomplete lineage sorting. Using our model, we propose methods for estimating the extent of the role of hybridization in both a likelihood and a Bayesian framework. The performance of our methods is examined using both simulated and empirical data.},
  Doi                      = {10.1016/j.tpb.2008.10.004},
  File                     = {:home/leo/articles/02-ImportedToMendeley/TheorPopBiol.2009.35doi10004.pdf:pdf},
  ISSN                     = {1096-0325},
  Keywords                 = {Animals,Bayes Theorem,Computer Simulation,Evolution,Gene Flow,Genetic,Genetic Speciation,Hybridization,Likelihood Functions,Models,Molecular,Phylogeny},
  Pmid                     = {19038278},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/19038278}
}

@Article{Yu2011,
  Title                    = {{Coalescent histories on phylogenetic networks and detection of hybridization despite incomplete lineage sorting.}},
  Author                   = {Yu, Yun and Than, Cuong and Degnan, James H and Nakhleh, Luay},
  Journal                  = {Systematic biology},
  Year                     = {2011},

  Month                    = mar,
  Number                   = {2},
  Pages                    = {138--49},
  Volume                   = {60},

  Abstract                 = {Analyses of the increasingly available genomic data continue to reveal the extent of hybridization and its role in the evolutionary diversification of various groups of species. We show, through extensive coalescent-based simulations of multilocus data sets on phylogenetic networks, how divergence times before and after hybridization events can result in incomplete lineage sorting with gene tree incongruence signatures identical to those exhibited by hybridization. Evolutionary analysis of such data under the assumption of a species tree model can miss all hybridization events, whereas analysis under the assumption of a species network model would grossly overestimate hybridization events. These issues necessitate a paradigm shift in evolutionary analysis under these scenarios, from a model that assumes a priori a single source of gene tree incongruence to one that integrates multiple sources in a unifying framework. We propose a framework of coalescence within the branches of a phylogenetic network and show how this framework can be used to detect hybridization despite incomplete lineage sorting. We apply the model to simulated data and show that the signature of hybridization can be revealed as long as the interval between the divergence times of the species involved in hybridization is not too small. We reanalyze a data set of 106 loci from 7 in-group Saccharomyces species for which a species tree with no hybridization has been reported in the literature. Our analysis supports the hypothesis that hybridization occurred during the evolution of this group, explaining a large amount of the incongruence in the data. Our findings show that an integrative approach to gene tree incongruence and its reconciliation is needed. Our framework will help in systematically analyzing genomic data for the occurrence of hybridization and elucidating its evolutionary role.},
  Doi                      = {10.1093/sysbio/syq084},
  File                     = {:home/leo/articles/02-ImportedToMendeley/SystBiol.2011.138syq084.pdf:pdf},
  ISSN                     = {1076-836X},
  Keywords                 = {Genes, Fungal,Hybridization, Genetic,Models, Genetic,Phylogeny,Saccharomyces,Saccharomyces: classification,Saccharomyces: genetics},
  Pmid                     = {21248369},
  Url                      = {http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3167682\&tool=pmcentrez\&rendertype=abstract}
}

@Article{Doyon2011,
  Title                    = {{An efficient algorithm for gene/species trees parsimonious reconciliation with losses, duplications and transfers}},
  Author                   = {Doyon, Jean-Philippe and Scornavacca, Celine and Gorbunov, K Yu and Sz{\"o}ll{\H{o}}si, Gergely J and Ranwez, Vincent and Berry, Vincent},
  Journal                  = {Comparative Genomics},
  Year                     = {2011},
  Number                   = {October},
  Pages                    = {93--108},

  File                     = {:home/leo/articles/02-ImportedToMendeley/LectureNotesCS.2010.093.pdf:PDF;:home/leo/articles/02-ImportedToMendeley/JPDoyon-RCG2010_slides.pdf:PDF;:home/leo/articles/02-ImportedToMendeley/JPDoyon-jobim2010.pdf:PDF;:home/leo/articles/02-ImportedToMendeley/JPDoyon-smbe2010.pdf:PDF;:home/leo/articles/02-ImportedToMendeley/Recomb.2010.93.pdf:PDF;:home/leo/articles/02-ImportedToMendeley/JPDoyon-smbe2010_report.pdf:PDF},
  Keywords                 = {DLT model,duplications,duploss,gene and species trees,losses,parsimony,reconciliation,transfers},
  Mendeley-tags            = {DLT model,duploss,reconciliation},
  Owner                    = {leo},
  Publisher                = {Springer},
  Timestamp                = {2013.10.14-1517},
  Url                      = {http://www.springerlink.com/index/D8H78633P3N31268.pdf}
}

@Article{Lanier2014,
  Title                    = {How low can you go? The effects of mutation rate on the accuracy of species-tree estimation},
  Author                   = {Lanier, Hayley C. and Huang, Huateng and Knowles, L. Lacey},
  Journal                  = {Molecular Phylogenetics and Evolution},
  Year                     = {2014},

  Month                    = {Sep},
  Pages                    = {112-119},
  Volume                   = {70},

  Doi                      = {10.1016/j.ympev.2013.09.006},
  File                     = {:home/leo/articles/021/MolPhylEvol.2013.doi09006.pdf:PDF},
  Owner                    = {leo},
  Publisher                = {Elsevier },
  Timestamp                = {2013.09.29},
  Url                      = {http://dx.doi.org/10.1016/j.ympev.2013.09.006}
}
@Article{Wu2012b,
  Title                    = {{Accounting For Alignment Uncertainty in Phylogenomics}},
  Author                   = {Wu, Martin and Chatterji, Sourav and Eisen, Jonathan a.},
  Journal                  = {PLoS ONE},
  Year                     = {2012},

  Month                    = jan,
  Number                   = {1},
  Pages                    = {e30288},
  Volume                   = {7},

  Doi                      = {10.1371/journal.pone.0030288},
  Editor                   = {Salemi, Marco},
  File                     = {:home/leo/articles/02-ImportedToMendeley/PLoSone.2012.0030288.pdf:pdf},
  ISSN                     = {1932-6203},
  Owner                    = {leo},
  Timestamp                = {2013.10.17-1823},
  Url                      = {http://dx.plos.org/10.1371/journal.pone.0030288}
}

@Article{Maddison1997,
  Title                    = {{Gene Trees in Species Trees}},
  Author                   = {Maddison, Wayne P},
  Journal                  = {Systematic Biology},
  Year                     = {1997},
  Number                   = {3},
  Pages                    = {523--536},
  Volume                   = {46},

  Doi                      = {10.1093/sysbio/46.3.523},
  File                     = {:home/leo/articles/02-ImportedToMendeley/SystBiol.1997.523.pdf:pdf;:home/leo/articles/02-ImportedToMendeley/rajan.pdf:pdf},
  ISSN                     = {1063-5157},
  Url                      = {http://sysbio.oxfordjournals.org/cgi/doi/10.1093/sysbio/46.3.523}
}

@Article{Bansal2010,
  Title                    = {{Robinson-Foulds supertrees.}},
  Author                   = {Bansal, Mukul S and Burleigh, J Gordon and Eulenstein, Oliver and Fern\'{a}ndez-Baca, David},
  Journal                  = {Algorithms for molecular biology},
  Year                     = {2010},

  Month                    = jan,
  Pages                    = {18},
  Volume                   = {5},

  Abstract                 = {Supertree methods synthesize collections of small phylogenetic trees with incomplete taxon overlap into comprehensive trees, or supertrees, that include all taxa found in the input trees. Supertree methods based on the well established Robinson-Foulds (RF) distance have the potential to build supertrees that retain much information from the input trees. Specifically, the RF supertree problem seeks a binary supertree that minimizes the sum of the RF distances from the supertree to the input trees. Thus, an RF supertree is a supertree that is consistent with the largest number of clusters (or clades) from the input trees.},
  Doi                      = {10.1186/1748-7188-5-18},
  File                     = {:home/leo/articles/02-ImportedToMendeley/AlgoMolBiol.2010.5-18.pdf:pdf},
  ISSN                     = {1748-7188},
  Pmid                     = {20181274},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/22431553}
}

@Article{Larget1999,
  Title                    = {{Markov chain Monte Carlo algorithms for the Bayesian analysis of phylogenetic trees}},
  Author                   = {Larget, B and Simon, DL},
  Journal                  = {Molecular Biology and Evolution},
  Year                     = {1999},
  Number                   = {6},
  Pages                    = {750--759},
  Volume                   = {16},

  File                     = {:home/leo/articles/02-ImportedToMendeley/MolBiolEvol.1999.750.pdf:pdf},
  Url                      = {http://www.renyi.hu/conferences/bp2008/Slides/larget1-budapest.pdf}
}
@Article{Yang1997,
  Title                    = {{Bayesian phylogenetic inference using DNA sequences: a Markov Chain Monte Carlo Method.}},
  Author                   = {Yang, Ziheng and Rannala, Bruce},
  Journal                  = {Molecular biology and evolution},
  Year                     = {1997},
  Number                   = {7},
  Pages                    = {717--24},
  Volume                   = {14},

  Abstract                 = {An improved Bayesian method is presented for estimating phylogenetic trees using DNA sequence data. The birth-death process with species sampling is used to specify the prior distribution of phylogenies and ancestral speciation times, and the posterior probabilities of phylogenies are used to estimate the maximum posterior probability (MAP) tree. Monte Carlo integration is used to integrate over the ancestral speciation times for particular trees. A Markov Chain Monte Carlo method is used to generate the set of trees with the highest posterior probabilities. Methods are described for an empirical Bayesian analysis, in which estimates of the speciation and extinction rates are used in calculating the posterior probabilities, and a hierarchical Bayesian analysis, in which these parameters are removed from the model by an additional integration. The Markov Chain Monte Carlo method avoids the requirement of our earlier method for calculating MAP trees to sum over all possible topologies (which limited the number of taxa in an analysis to about five). The methods are applied to analyze DNA sequences for nine species of primates, and the MAP tree, which is identical to a maximum-likelihood estimate of topology, has a probability of approximately 95\%.},
  File                     = {:home/leo/articles/02-ImportedToMendeley/MolBiolEvol.1997.717.pdf:pdf},
  ISSN                     = {0737-4038},
  Keywords                 = {Algorithms,Animals,Bayes Theorem,DNA,Mitochondrial,Mitochondrial: genetics,Monte Carlo Method,Phylogeny,Primates,Primates: genetics,Sequence Analysis,Sequence Analysis: methods},
  Pmid                     = {9214744},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/9214744}
}

@Article{Guigo1996,
  Title                    = {{Reconstruction of ancient molecular phylogeny.}},
  Author                   = {Guig\'{o}, R and Muchnik, I and Smith, T F},
  Journal                  = {Molecular phylogenetics and evolution},
  Year                     = {1996},

  Month                    = oct,
  Number                   = {2},
  Pages                    = {189--213},
  Volume                   = {6},

  Abstract                 = {Support for contradictory phylogenies is often obtained when molecular sequence data from different genes is used to reconstruct phylogenetic histories. Contradictory phylogenies can result from many data anomalies including unrecognized paralogy. Paralogy, defined as the reconstruction of a phylogenetic tree from a mixture of genes generated by duplications, has generally not been formally included in phylogenetic reconstructions. Here we undertake the task of reconstructing a single most likely evolutionary relationship among a range of taxa from a large set of apparently inconsistent gene trees. Under the assumption that differences among gene trees can be explained by gene duplications, and consequent losses, we have developed a method to obtain the global phylogeny minimizing the total number of postulated duplications and losses and to trace back such individual gene duplications to global genome duplications. We have used this method to infer the most likely phylogenetic relationship among 16 major higher eukaryotic taxa from the sequences of 53 different genes. Only five independent genome duplication events need to be postulated in order to explain the inconsistencies among these trees.},
  Doi                      = {10.1006/mpev.1996.0071},
  File                     = {:home/leo/articles/02-ImportedToMendeley/MolPhylogEvol.1996.189.pdf:pdf},
  ISSN                     = {1055-7903},
  Keywords                 = {Algorithms,Animals,Biological Evolution,Genes,Models, Biological,Multigene Family,Phylogeny,Species Specificity,duploss,phylogenomics,reconciliation},
  Mendeley-tags            = {duploss,phylogenomics,reconciliation},
  Pmid                     = {8899723},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/8899723}
}

@Article{Zhang1997,
  Title                    = {{On a Mirkin-Muchnik-Smith conjecture for comparing molecular phylogenies}},
  Author                   = {Zhang, L},
  Journal                  = {Journal of Computational Biology},
  Year                     = {1997},
  Number                   = {2},
  Pages                    = {177--187},
  Volume                   = {4},

  File                     = {:home/leo/articles/02-ImportedToMendeley/JComputBiol.1997.177.pdf:pdf},
  Url                      = {http://online.liebertpub.com/doi/abs/10.1089/cmb.1997.4.177}
}


@Article{Page1994,
  Title                    = {maps between trees and cladistic analysis of historical associations among genes, organisms and areas},
  Author                   = {Page, RDM},
  Journal                  = {Systematic Biology},
  Year                     = {1994},
  Number                   = {1},
  Pages                    = {58--77},
  Volume                   = {43},

  File                     = {:home/leo/articles/02-ImportedToMendeley/SystBiol.1994.58.pdf:pdf},
  Url                      = {http://sysbio.oxfordjournals.org/content/43/1/58.short}
}

%% review of GTP, talk about trimmed vc untrimmed species trees (ours are trimmed)
@Article{Bansal13, 
  Title                    = {Algorithms for Genome-Scale Phylogenetics Using Gene Tree Parsimony},
  Author                   = {Mukul S. Bansal and Oliver Eulenstein},
  Year                     = {2013},
  journal={Computational Biology and Bioinformatics, IEEE/ACM Transactions on},
  volume={10},
  number={4},
  pages={939--956},
  publisher={IEEE},
  Abstract                 = {The use of genomic datasets for phylogenetics is complicated by the fact that evolutionary processes such as gene duplication and loss, or incomplete lineage sorting (deep coalescence) cause incongruence among gene trees. One well- known approach that deals with this complication is gene tree parsimony, which, given a collection of gene trees, seeks a species tree that requires the smallest number of evolutionary events to explain the incongruence of the gene trees. However, a lack of efficient algorithms has limited the use of this approach. Here, we present efficient algorithms for SPR and TBR based local search heuristics for gene tree parsimony under the (i) duplication, (ii) loss, (iii) duplication-loss, and (iv) deep coalescence reconciliation costs. These novel algorithms improve upon the time complexities of previous algorithms for these problems by a factor of n, where n is the number of species in the collection of gene trees. Our algorithms provide a substantial improvement in runtime and scalability compared to previous implementations and enable large-scale gene tree parsimony analyses using any of the four reconciliation costs. Our algorithms have been implemented in the software packages DupTree and iGTP, and have already been used to perform several compelling phylogenetic studies.},
  File                     = {:home/leo/articles/021/Bansal_TCBB2013.pdf:PDF},
  Owner                    = {leo},
  Timestamp                = {2013.09.29}
}

@Article{Wu2011,
  Title                    = {{Structural properties of the reconciliation space and their applications in enumerating nearly-optimal reconciliations between a gene tree and a species tree.}},
  Author                   = {Wu, Taoyang and Zhang, Louxin},
  Journal                  = {BMC bioinformatics},
  Year                     = {2011},

  Month                    = jan,
  Number                   = {Suppl 9},
  Pages                    = {S7},
  Volume                   = {12 Suppl 9},

  Abstract                 = {A gene tree for a gene family is often discordant with the containing species tree because of its complex evolutionary course during which gene duplication, gene loss and incomplete lineage sorting events might occur. Hence, it is of great challenge to infer the containing species tree from a set of gene trees. One common approach to this inference problem is through gene tree and species tree reconciliation.},
  Doi                      = {10.1186/1471-2105-12-S9-S7},
  File                     = {:home/leo/articles/02-ImportedToMendeley/BMCBioinfo.2011.12-S9-S7.pdf:pdf},
  ISSN                     = {1471-2105},
  Keywords                 = {duploss,phylogenomics,reconciliation},
  Mendeley-tags            = {duploss,phylogenomics,reconciliation},
  Pmid                     = {22151151},
  Publisher                = {BioMed Central Ltd},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/22151151}
}

%original RF
@Article{Robinson1981,
  Title                    = {{Comparison of phylogenetic trees}},
  Author                   = {Robinson, DF and Foulds, LR},
  Journal                  = {Mathematical Biosciences},
  Year                     = {1981},
  Pages                    = {131--147},
  Volume                   = {53},

  File                     = {:home/leo/articles/02-ImportedToMendeley/MathBiosci.1981.131.pdf:pdf},
  Url                      = {http://www.sciencedirect.com/science/article/pii/0025556481900432}
}
% rooted RF
@Article{Day1985,
  Title                    = {{Optimal algorithms for comparing trees with labeled leaves}},
  Author                   = {Day, WHE},
  Journal                  = {Journal of Classification},
  Year                     = {1985},
  Pages                    = {7--28},
  Volume                   = {2},

  File                     =  {:home/leo/articles/02-ImportedToMendeley/JournalClassification.1985.7.pdf:pdf},
  Keywords                 = {algorithm complexity,algorithm design,classifications,comparing hierarchical,comparing phylogenetic trees,consensus index,strict con-},
  Url                      = {http://www.springerlink.com/index/q5906x80g44p44k8.pdf}
}
% rooted with unrooted RF
@Article{Pawe2012,
  Title                    = {{A Robinson-Foulds Measure to Compare Unrooted Trees with Rooted Trees}},
  Author                   = {G\'{o}recki, Pawe and Eulenstein, Oliver},
  Journal                  = {Bioinformatics Research and Applications},
  Year                     = {2012},
  Pages                    = {115--126},
  Volume                   = {7292},

  Abstract                 = {The Robinson-Foulds (RF) distance is a well-established measure for the comparative analysis of phylogenetic trees, and comparing unrooted gene trees with rooted species trees is a standard application in phylogenetics. However, the RF distance is not defined to compare unrooted trees with rooted trees. Here we propose urRF a new measure based on the standard RF distance to compare an unrooted tree with a rooted tree, and describe a linear time algorithm to compute urRF. Further, we show that our measure relates to several intriguing properties that can be crucial for the in-depth comparative study and the credible rooting of trees. Finally, we present detailed empirical studies based on real data sets demonstrating the performance of our novel RF measure, the efficiency of our linear time algorithm and the quality of determining optimal rootings of unrooted gene trees.},
  Doi                      = {10.1007/978-3-642-30191-9\_12},
  File                     = {:home/leo/articles/02-ImportedToMendeley/LectureNotesCS.2012.115.pdf:pdf},
  Keywords                 = {00.distance},
  Owner                    = {leo},
  Timestamp                = {20131130.2357},
  Url                      = {http://link.springer.com/chapter/10.1007\%2F978-3-642-30191-9\_12}
}

% Hdist
@Article{Nye2006,
  Title                    = {{A novel algorithm and web-based tool for comparing two alternative phylogenetic trees.}},
  Author                   = {Nye, Tom M W and Li\`{o}, Pietro and Gilks, Walter R},
  Journal                  = {Bioinformatics},
  Year                     = {2006},

  Month                    = jan,
  Number                   = {1},
  Pages                    = {117--9},
  Volume                   = {22},

  Abstract                 = {SUMMARY: We describe an algorithm and software tool for comparing alternative phylogenetic trees. The main application of the software is to compare phylogenies obtained using different phylogenetic methods for some fixed set of species or obtained using different gene sequences from those species. The algorithm pairs up each branch in one phylogeny with a matching branch in the second phylogeny and finds the optimum 1-to-1 map between branches in the two trees in terms of a topological score. The software enables the user to explore the corresponding mapping between the phylogenies interactively, and clearly highlights those parts of the trees that differ, both in terms of topology and branch length. AVAILABILITY: The software is implemented as a Java applet at http://www.mrc-bsu.cam.ac.uk/personal/thomas/phylo\_comparison/comparison\_page.html. It is also available on request from the authors.},
  Doi                      = {10.1093/bioinformatics/bti720},
  File                     = {:home/leo/articles/02-ImportedToMendeley/Bioinformatics.2006.117.pdf:pdf},
  ISSN                     = {1367-4803},
  Keywords                 = {Algorithms,Computational Biology,Computational Biology: methods,Computer Graphics,Genetic,HIV,HIV: genetics,Internet,Models,Phylogeny,Programming Languages,Sequence Alignment,Software,Statistical,User-Computer Interface},
  Pmid                     = {16234319},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/16234319}
}
% extention of Nye2006 metric
@Article{Bogdanowicz08,
  Title                    = {Comparing phylogenetic trees using a minimum weight perfect matching},
  Author                   = {Bogdanowicz, Damian},
  Pages                    = {1-4},
  Year                     = {2008},
  Month                    = {May},

  Journal                = {2008 1st International Conference on Information Technology},
  Doi                      = {10.1109/INFTECH.2008.4621680},
  File                     = {:home/leo/articles/021/ProcIntlConfInfoTech.2008.1.pdf:PDF},
  ISBN                     = {978-1-4244-2244-9},
  Keywords                 = {00.distance},
  Owner                    = {leo},
  Timestamp                = {20140113.1823},
  Url                      = {http://dx.doi.org/10.1109/INFTECH.2008.4621680}
}

% SPR distance
@Article{Beiko2006,
  Title                    = {{Phylogenetic identification of lateral genetic transfer events.}},
  Author                   = {Beiko, Robert G and Hamilton, Nicholas},
  Journal                  = {BMC Evolutionary Biology},
  Year                     = {2006},

  Month                    = jan,
  Pages                    = {15},
  Volume                   = {6},

  Abstract                 = {Lateral genetic transfer can lead to disagreements among phylogenetic trees comprising sequences from the same set of taxa. Where topological discordance is thought to have arisen through genetic transfer events, tree comparisons can be used to identify the lineages that may have shared genetic information. An 'edit path' of one or more transfer events can be represented with a series of subtree prune and regraft (SPR) operations, but finding the optimal such set of operations is NP-hard for comparisons between rooted trees, and may be so for unrooted trees as well.},
  Doi                      = {10.1186/1471-2148-6-15},
  File                     = {:home/leo/articles/02-ImportedToMendeley/BMCEvolBiol.6-15.pdf:pdf},
  ISSN                     = {1471-2148},
  Keywords                 = {Algorithms,Computational Biology,Computational Biology: methods,Evolution, Molecular,Gene Transfer, Horizontal,Genome, Archaeal,Genome, Bacterial,Genomics,Genomics: methods,Phylogeny},
  Pmid                     = {16472400},
  Url                      = {http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=1431587\&tool=pmcentrez\&rendertype=abstract}
}


@book{Liu2001,
  title={Monte Carlo strategies in scientific computing},
  author={Liu, Jun S},
  year={2001},
  publisher={Springer}
}

% matrix repres parsimony (MRP) 
@Article{Bininda2004,
  Title                    = {{The evolution of supertrees.}},
  Author                   = {Bininda-Emonds, Olaf R P},
  Journal                  = {Trends in Ecology \& Evolution},
  Year                     = {2004},
  Month                    = jun,
  Number                   = {6},
  Pages                    = {315--22},
  Volume                   = {19},
  Abstract                 = {Supertrees result from combining many smaller, overlapping phylogenetic trees into a single, more comprehensive tree. As such, supertree construction is probably as old as the field of systematics itself, and remains our only way of visualizing the Tree of Life as a whole. Over the past decade, supertree construction has gained a more formal, objective footing, and has become an area of active theoretical and practical research. Here, I review the history of the supertree approach, focusing mainly on its current implementation. The supertrees of today represent some of the largest, complete phylogenies available for many groups, but are not without their critics. I conclude by arguing that the ever-growing molecular revolution will result in supertree construction taking on a new role and implementation in the future for analyzing large DNA sequence matrices as part of a divide-and-conquer phylogenetic approach.},
  Doi                      = {10.1016/j.tree.2004.03.015},
  File                     = {:home/leo/articles/02-ImportedToMendeley/TrendsEcolEvol.2004.315.pdf:pdf},
  ISSN                     = {0169-5347},
  Pmid                     = {16701277},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/16701277}
}
% MRP
@Article{Ragan1992,
  Title                    = {Phylogenetic inference based on matrix representation of trees.},
  Author                   = {Ragan, M. A.},
  Journal                  = {Mol Phylogenet Evol},
  Year                     = {1992},
  Month                    = {Mar},
  Number                   = {1},
  Pages                    = {53--58},
  Volume                   = {1},
  Abstract                 = {Rooted phylogenetic trees can be represented as matrices in which the rows correspond to termini, and columns correspond to internal nodes (elements of the n-tree). Parsimony analysis of such a matrix will fully recover the topology of the original tree. The maximum size of the represented matrix depends only on the number of termini in the tree; for a tree derived from molecular sequences, the represented matrix may be orders of magnitude smaller than the original data matrix. Representations of multiple trees (which may or may not have identical termini) can readily be combined into a single matrix; columns of discrete-character-state data can be added and, if desired, weighted differentially. Parsimony analysis of the resulting composite matrix yields a hybrid supertree which typically provides greater resolution than conventional consensus trees. Use of this method is illustrated with examples involving multiple tRNA genes in organelles and multiple protein-coding genes in eukaryotes.},
  File                     = {:home/leo/articles/021/MolPhylogEvol.1992.53.pdf:PDF},
  Institution              = {Institute for Marine Biosciences, National Research Council of Canada, Halifax, Nova Scotia.},
  Keywords                 = {Animals; DNA, genetics; Dictyostelium, classification/genetics; Genetic Techniques; Humans; Microcomputers; Organelles, metabolism; Phylogeny; RNA, Transfer, genetics},
  Medline-pst              = {ppublish},
  Owner                    = {leo},
  Pmid                     = {1342924},
  Timestamp                = {20140113.1446},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/1342924}
}

@Article{smith1992,
  Title                    = {Bayesian statistics without tears: a sampling--resampling perspective},
  Author                   = {Smith, Adrian FM and Gelfand, Alan E},
  Journal                  = {The American Statistician},
  Year                     = {1992},
  Number                   = {2},
  Pages                    = {84--88},
  Volume                   = {46},

  File                     = {:home/leo/articles/021/AmStats.1992.84.pdf:PDF;:home/leo/articles/021/AmStats.1992.84_notes.pdf:PDF},
  Keywords                 = {00.importancesampling},
  Publisher                = {Taylor \& Francis Group},
  Timestamp                = {20140107.1615},
  Url                      = {http://www.jstor.org/stable/2684170}
}

@Article{Chung2013,
  Title                    = {Computing the Joint Distribution of Tree Shape and Tree Distance for Gene Tree Inference and Recombination Detection},
  Author                   = {Chung, Yujin and Perna, Nicole T. and Ane, Cecile},
  Journal                  = {IEEE/ACM Transactions on Computational Biology and Bioinformatics},
  Year                     = {2013},
 Note = {preprint},
  Doi                      = {10.1109/TCBB.2013.109},
  File                     = {:home/leo/articles/021/TCBB.2013.109.pdf:PDF},
  Owner                    = {leo},
  Publisher                = {Institute of Electrical and Electronics Engineers},
  Timestamp                = {20140114.0303},
  Url                      = {http://dx.doi.org/10.1109/TCBB.2013.109}
}

% pseudorandom parallel
@Article{Feng2003,
  Title                    = {Parallel algorithms for {Bayes}ian phylogenetic inference},
  Author                   = {Feng, Xizhou and Buell, Duncan A. and Rose, John R. and Waddell, Peter J.},
  Journal                  = {Journal of Parallel and Distributed Computing},
  Year                     = {2003},
 File                     = {:home/leo/articles/021/JParallelDistrComp.2003.707.pdf:PDF},
  Month                    = {Jul},
  Number                   = {7-8},
  Pages                    = {707-718},
  Volume                   = {63},

  Doi                      = {10.1016/S0743-7315(03)00079-0},
  Owner                    = {leo},
  Publisher                = {Elsevier },
  Timestamp                = {20140114.1830},
  Url                      = {http://dx.doi.org/10.1016/S0743-7315(03)00079-0}
}
@Article{Feng2006,
  Title                    = {PBPI: a high performance implementation of Bayesian phylogenetic inference},
  Author                   = {Feng, Xizhou and Cameron, Kirk W. and Buell, Duncan A.},
  Note                    = {Article no 75},
  Publisher                = {Institute of Electrical and Electronics Engineers},
  Year                     = {2006},

  Journal                = {Proceedings of the 2006 ACM/IEEE conference on Supercomputing},
  Doi                      = {10.1145/1188455.1188535},
  File                     = {:home/leo/articles/021/ProcACMConferenceSupercomput.2006.292.pdf:PDF},
  ISBN                     = {0769527000},
  Owner                    = {leo},
  Timestamp                = {20140114.1832},
  Url                      = {http://dx.doi.org/10.1145/1188455.1188535}
}

@article{Martins08,
  pmid = {18612422},
  author = {de Oliveira Martins, Leonardo and Leal, Elcio and Kishino, Hirohisa},
  title = {Phylogenetic detection of recombination with a Bayesian prior on the distance between trees.},
  journal = {PLoS ONE},
  volume = {3},
  number = {7},
  pages = {e2651},
  year = {2008},
  doi = {10.1371/journal.pone.0002651},
  file = {file:///home/leo/articles/personal.publications/PLoSone.0002651.pdf:pdf},
}
@article{Martins10,
  Publisher = {Springer},
  author = {de Oliveira Martins, Leonardo and Kishino, Hirohisa},
  title = {Distribution of distances between topologies and its effect on detection of phylogenetic recombination},
  journal = {Annals of the Institute of Statistical Mathematics},
  volume = {62},
  number = {1},
  pages = {145-159},
  year = {2010},
  doi = {10.1007/s10463-009-0259-8},
  file = {file:///home/leo/articles/personal.publications/AnnInstStatMaths.2010.145.pdf:pdf},
}

@Article{Lartillot2004,
  Title                    = {{A Bayesian mixture model for across-site heterogeneities in the amino-acid replacement process.}},
  Author                   = {Lartillot, Nicolas and Philippe, Herv\'{e}},
  Journal                  = {Molecular biology and evolution},
  Year                     = {2004},

  Month                    = jun,
  Number                   = {6},
  Pages                    = {1095--109},
  Volume                   = {21},

  Abstract                 = {Most current models of sequence evolution assume that all sites of a protein evolve under the same substitution process, characterized by a 20 x 20 substitution matrix. Here, we propose to relax this assumption by developing a Bayesian mixture model that allows the amino-acid replacement pattern at different sites of a protein alignment to be described by distinct substitution processes. Our model, named CAT, assumes the existence of distinct processes (or classes) differing by their equilibrium frequencies over the 20 residues. Through the use of a Dirichlet process prior, the total number of classes and their respective amino-acid profiles, as well as the affiliations of each site to a given class, are all free variables of the model. In this way, the CAT model is able to adapt to the complexity actually present in the data, and it yields an estimate of the substitutional heterogeneity through the posterior mean number of classes. We show that a significant level of heterogeneity is present in the substitution patterns of proteins, and that the standard one-matrix model fails to account for this heterogeneity. By evaluating the Bayes factor, we demonstrate that the standard model is outperformed by CAT on all of the data sets which we analyzed. Altogether, these results suggest that the complexity of the pattern of substitution of real sequences is better captured by the CAT model, offering the possibility of studying its impact on phylogenetic reconstruction and its connections with structure-function determinants.},
  Doi                      = {10.1093/molbev/msh112},
  File                     = {:home/leo/articles/02-ImportedToMendeley/MolBiolEvol.2004.1095.pdf:pdf},
  ISSN                     = {0737-4038},
  Keywords                 = {Amino Acid Sequence,Amino Acid Substitution,Amino Acid Substitution: genetics,Bayes Theorem,Evolution, Molecular,Models, Genetic,Phylogeny,Sequence Alignment},
  Owner                    = {leo},
  Pmid                     = {15014145},
  Timestamp                = {2013.10.14-1609},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/15014145}
}
% simmul anneal
@Article{Rubenthaler09,
  Title                    = {Fast simulated annealing in Rd with an application to maximum likelihood estimation in state-space models},
  Author                   = {Rubenthaler, Sylvain and Ryd{\'e}n, Tobias and Wiktorsson, Magnus},
  Journal                  = {Stochastic Processes and their Applications},
  Year                     = {2009},
  Number                   = {6},
  Pages                    = {1912--1931},
  Volume                   = {119},

  Arxivid                  = {math.PR/0609353},
  Doi                      = {10.1016/j.spa.2008.09.007},
  Eprint                   = {0609353},
  File                     = {:home/leo/articles/021/StochProcessesAppl.2009.1912.pdf:PDF},
  Publisher                = {Elsevier},
  Timestamp                = {20140114.2200}
}

@Article{Fletcher2009,
  Title                    = {{INDELible: a flexible simulator of biological sequence evolution.}},
  Author                   = {Fletcher, William and Yang, Ziheng},
  Journal                  = {Molecular biology and evolution},
  Year                     = {2009},

  Month                    = aug,
  Number                   = {8},
  Pages                    = {1879--88},
  Volume                   = {26},

  Abstract                 = {Many methods exist for reconstructing phylogenies from molecular sequence data, but few phylogenies are known and can be used to check their efficacy. Simulation remains the most important approach to testing the accuracy and robustness of phylogenetic inference methods. However, current simulation programs are limited, especially concerning realistic models for simulating insertions and deletions. We implement a portable and flexible application, named INDELible, for generating nucleotide, amino acid and codon sequence data by simulating insertions and deletions (indels) as well as substitutions. Indels are simulated under several models of indel-length distribution. The program implements a rich repertoire of substitution models, including the general unrestricted model and nonstationary nonhomogeneous models of nucleotide substitution, mixture, and partition models that account for heterogeneity among sites, and codon models that allow the nonsynonymous/synonymous substitution rate ratio to vary among sites and branches. With its many unique features, INDELible should be useful for evaluating the performance of many inference methods, including those for multiple sequence alignment, phylogenetic tree inference, and ancestral sequence, or genome reconstruction.},
  Doi                      = {10.1093/molbev/msp098},
  File                     = {:home/leo/articles/02-ImportedToMendeley/MolBiolEvol.2009.1879.pdf:pdf},
  ISSN                     = {1537-1719},
  Keywords                 = {Algorithms,Animals,Evolution, Molecular,Humans,INDEL Mutation,Models, Genetic,Software},
  Owner                    = {leo},
  Pmid                     = {19423664},
  Timestamp                = {2013.05.04},
  Url                      = {http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2712615\&tool=pmcentrez\&rendertype=abstract}
}

@Article{Yang2007,
  Title                    = {{PAML 4: phylogenetic analysis by maximum likelihood.}},
  Author                   = {Yang, Ziheng},
  Journal                  = {Molecular biology and evolution},
  Year                     = {2007},

  Month                    = aug,
  Number                   = {8},
  Pages                    = {1586--91},
  Volume                   = {24},

  Abstract                 = {PAML, currently in version 4, is a package of programs for phylogenetic analyses of DNA and protein sequences using maximum likelihood (ML). The programs may be used to compare and test phylogenetic trees, but their main strengths lie in the rich repertoire of evolutionary models implemented, which can be used to estimate parameters in models of sequence evolution and to test interesting biological hypotheses. Uses of the programs include estimation of synonymous and nonsynonymous rates (d(N) and d(S)) between two protein-coding DNA sequences, inference of positive Darwinian selection through phylogenetic comparison of protein-coding genes, reconstruction of ancestral genes and proteins for molecular restoration studies of extinct life forms, combined analysis of heterogeneous data sets from multiple gene loci, and estimation of species divergence times incorporating uncertainties in fossil calibrations. This note discusses some of the major applications of the package, which includes example data sets to demonstrate their use. The package is written in ANSI C, and runs under Windows, Mac OSX, and UNIX systems. It is available at -- (http://abacus.gene.ucl.ac.uk/software/paml.html).},
  Doi                      = {10.1093/molbev/msm088},
  File                     = {:home/leo/articles/02-ImportedToMendeley/MolBiolEvol.2007.1586.pdf:pdf},
  ISSN                     = {0737-4038},
  Keywords                 = {Animals,Computer Simulation,Genetic Variation,Likelihood Functions,Models, Genetic,Phylogeny,Selection, Genetic,Software,Species Specificity},
  Owner                    = {leo},
  Pmid                     = {17483113},
  Timestamp                = {2013.05.04},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/17483113}
}
@Article{Drummond2012,
  Title                    = {{Bayesian phylogenetics with BEAUti and the BEAST 1.7.}},
  Author                   = {Drummond, Alexei J. and Suchard, Marc A. and Xie, Dong and Rambaut, Andrew},
  Journal                  = {Molecular Biology and Evolution},
  Year                     = {2012},

  Month                    = {Aug},
  Number                   = {8},
  Pages                    = {1969--1973},
  Volume                   = {29},

  Abstract                 = {Computational evolutionary biology, statistical phylogenetics and coalescent-based population genetics are becoming increasingly central to the analysis and understanding of molecular sequence data. We present the Bayesian Evolutionary Analysis by Sampling Trees (BEAST) software package version 1.7, which implements a family of Markov chain Monte Carlo (MCMC) algorithms for Bayesian phylogenetic inference, divergence time dating, coalescent analysis, phylogeography and related molecular evolutionary analyses. This package includes an enhanced graphical user interface program called Bayesian Evolutionary Analysis Utility (BEAUti) that enables access to advanced models for molecular sequence and phenotypic trait evolution that were previously available to developers only. The package also provides new tools for visualizing and summarizing multispecies coalescent and phylogeographic analyses. BEAUti and BEAST 1.7 are open source under the GNU lesser general public license and available at http://beast-mcmc.googlecode.com and http://beast.bio.ed.ac.uk.},
  Doi                      = {10.1093/molbev/mss075},
  Institution              = {Department of Computer Science, University of Auckland, Auckland, New Zealand. alexei@cs.auckland.ac.nz},
  Keywords                 = {Animals; Base Sequence; Bayes Theorem; Computational Biology, methods; DNA, Mitochondrial, genetics; Finches, genetics; Molecular Sequence Data; Phenotype; Phylogeny; Software; User-Computer Interface},
  Medline-pst              = {ppublish},
  Owner                    = {leo},
  Pii                      = {mss075},
 File                     = {:home/leo/articles/02-ImportedToMendeley/MolBiolEvol.2012.mss075.pdf:pdf},
  Pmid                     = {22367748},
  Timestamp                = {20140205.0038},
 ISSN                     = {1537-1719},
  Url                      = {http://dx.doi.org/10.1093/molbev/mss075}
}
% NNI
@Article{Allen2001,
  Title                    = {{Subtree Transfer Operations and Their Induced Metrics on Evolutionary Trees}},
  Author                   = {Allen, Benjamin L. and Steel, Mike},
  Journal                  = {Annals of Combinatorics},
  Year                     = {2001},

  Month                    = jun,
  Number                   = {1},
  Pages                    = {1--15},
  Volume                   = {5},

  Doi                      = {10.1007/s00026-001-8006-8},
  File                     = {:home/leo/articles/02-ImportedToMendeley/AnnalsCombinatorics.2001.1.pdf:pdf},
  ISSN                     = {0218-0006},
  Keywords                 = {fixed parameter tractability,metrics,subtree transfer,trees},
  Url                      = {http://www.springerlink.com/openurl.asp?genre=article\&id=doi:10.1007/s00026-001-8006-8}
}
@Article{Waterman1978,
  Title                    = {On the similarity of dendrograms},
  Author                   = {Waterman, M.S. and Smith, T.F.},
  Journal                  = {Journal of Theoretical Biology},
  Year                     = {1978},
 File                     = {:home/leo/articles/021/JTheorBiol.1978.789.pdf:PDF},
  Month                    = {Aug},
  Number                   = {4},
  Pages                    = {789-800},
  Volume                   = {73},

  Doi                      = {10.1016/0022-5193(78)90137-6},
  Owner                    = {leo},
  Publisher                = {Elsevier },
  Timestamp                = {20140115.1745},
  Url                      = {http://dx.doi.org/10.1016/0022-5193(78)90137-6}
}
@Article{Maddison2006,
  Title                    = {{Inferring phylogeny despite incomplete lineage sorting.}},
  Author                   = {Maddison, Wayne P and Knowles, L Lacey},
  Journal                  = {Systematic biology},
  Year                     = {2006},

  Month                    = feb,
  Number                   = {1},
  Pages                    = {21--30},
  Volume                   = {55},

  Abstract                 = {It is now well known that incomplete lineage sorting can cause serious difficulties for phylogenetic inference, but little attention has been paid to methods that attempt to overcome these difficulties by explicitly considering the processes that produce them. Here we explore approaches to phylogenetic inference designed to consider retention and sorting of ancestral polymorphism. We examine how the reconstructability of a species (or population) phylogeny is affected by (a) the number of loci used to estimate the phylogeny and (b) the number of individuals sampled per species. Even in difficult cases with considerable incomplete lineage sorting (times between divergences less than 1 N(e) generations), we found the reconstructed species trees matched the "true" species trees in at least three out of five partitions, as long as a reasonable number of individuals per species were sampled. We also studied the tradeoff between sampling more loci versus more individuals. Although increasing the number of loci gives more accurate trees for a given sampling effort with deeper species trees (e.g., total depth of 10 N(e) generations), sampling more individuals often gives better results than sampling more loci with shallower species trees (e.g., depth = 1 N(e)). Taken together, these results demonstrate that gene sequences retain enough signal to achieve an accurate estimate of phylogeny despite widespread incomplete lineage sorting. Continued improvement in our methods to reconstruct phylogeny near the species level will require a shift to a compound model that considers not only nucleotide or character state substitutions, but also the population genetics processes of lineage sorting. [Coalescence; divergence; population; speciation.].},
  Doi                      = {10.1080/10635150500354928},
  File                     = {:home/leo/articles/02-ImportedToMendeley/SystBiol.2006.21.pdf:pdf},
  ISSN                     = {1063-5157},
  Keywords                 = {Evolution, Molecular,Genetic Speciation,Genetics, Population,Phylogeny,Polymorphism, Genetic,Selection Bias},
  Pmid                     = {16507521},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/16507521}
}

@Article{Steel1993,
  Title                    = {{Distributions of Tree Comparison Metrics--Some New Results}},
  Author                   = {Steel, M. A. and Penny, David},
  Journal                  = {Systematic Biology},
  Year                     = {1993},

  Month                    = jun,
  Number                   = {2},
  Pages                    = {126--141},
  Volume                   = {42},

  Doi                      = {10.1093/sysbio/42.2.126},
  File                     = {:home/leo/articles/02-ImportedToMendeley/SystBiol.1993.126.pdf:pdf},
  ISSN                     = {1063-5157},
  Url                      = {http://sysbio.oxfordjournals.org/cgi/doi/10.1093/sysbio/42.2.126}
}

@Article{Wu2013,
  Title                    = {Most parsimonious reconciliation in the presence of gene duplication, loss, and deep coalescence using labeled coalescent trees},
  Author                   = {Wu, Y.-C. and Rasmussen, M. D. and Bansal, M. S. and Kellis, M.},
  Journal                  = {Genome Research},
  Year                     = {2013},
  note = {doi: 10.1101/gr.161968.113},
  Month                    = {Dec},

  Doi                      = {10.1101/gr.161968.113},
  File                     = {:home/leo/articles/021/GenomeRes.2013.gr161968113.pdf:PDF},
  Keywords                 = {00.reconciliation},
  Owner                    = {leo},
  Publisher                = {Cold Spring Harbor Laboratory Press},
  
  Timestamp                = {20131208.1430},
  Url                      = {http://dx.doi.org/10.1101/gr.161968.113}
}

@Article{Stark2007,
  Title                    = {Discovery of functional elements in 12 Drosophila genomes using evolutionary signatures},
  Author                   = {Stark, Alexander and Lin, Michael F. and Kheradpour, Pouya and Pedersen, Jakob S. and Parts, Leopold and Carlson, Joseph W. and Crosby, Madeline A. and Rasmussen, Matthew D. and Roy, Sushmita and Deoras, Ameya N. and et al.},
  Journal                  = {Nature},
  Year                     = {2007},

  Month                    = {Nov},
  Number                   = {7167},
  Pages                    = {219-232},
  Volume                   = {450},

  Doi                      = {10.1038/nature06340},
  File                     = {:home/leo/articles/021/Nature.2007.219.pdf:PDF},
  Keywords                 = {00.drosophila},
  Owner                    = {leo},
  Publisher                = {Nature Publishing Group},
  Timestamp                = {20131107.1555},
  Url                      = {http://dx.doi.org/10.1038/nature06340}
}
## duplication-aware phylogenomic analyses
@Article{Katz2012,
  Title                    = {Turning the Crown Upside Down: Gene Tree Parsimony Roots the Eukaryotic Tree of Life},
  Author                   = {Katz, L. A. and Grant, J. R. and Parfrey, L. W. and Burleigh, J. G.},
  Journal                  = {Systematic Biology},
  Year                     = {2012},

  Month                    = {Jun},
  Number                   = {4},
  Pages                    = {653-660},
  Volume                   = {61},

  Doi                      = {10.1093/sysbio/sys026},
  File                     = {:home/leo/articles/02-ImportedToMendeley/SystBiol.2012.sys026.pdf:pdf},
  ISSN                     = {1063-5157},
  Keywords                 = {00.reconciliation},
  Owner                    = {leo},
  Publisher                = {Oxford University Press},
  Timestamp                = {20140205.0047},
  Url                      = {http://dx.doi.org/10.1093/sysbio/sys026}
}

## duplication-aware phylogenomic analyses
@Article{Holton2010,
  Title                    = {{Deep genomic-scale analyses of the metazoa reject Coelomata: evidence from single- and multigene families analyzed under a supertree and supermatrix paradigm.}},
  Author                   = {Holton, Th\'{e}r\`{e}se a and Pisani, Davide},
  Journal                  = {Genome biology and evolution},
  Year                     = {2010},

  Month                    = jan,
  Pages                    = {310--24},
  Volume                   = {2},

  Abstract                 = {Solving the phylogeny of the animals with bilateral symmetry has proven difficult. Morphological studies have suggested a variety of alternative hypotheses, of which, Hyman's Coelomata hypothesis has become the most established. Studies based on 18S rRNA have failed to endorse Coelomata, supporting instead the rearrangement of the protostomes into two new clades: the Lophotrochozoa (including, e.g., the molluscs and the annelids) and the Ecdysozoa (including the Panarthropoda and most pseudocoelomates, such as the nematodes and priapulids). Support for this new animal phylogeny has been attained from expressed sequence tag studies, although these generally have a limited gene sampling. In contrast, deep genomic-scale analyses have often supported Coelomata. However, these studies are problematic due to their limited taxonomic sampling, which could exacerbate tree reconstruction artifacts. Here, we address both of these sampling limitations; we study the effect of long-branch attraction (LBA) in deep genomic-scale analyses and provide convincing evidence, using both single- and multigene families, that Coelomata is an artifact. We show that optimal outgroup selection is key in avoiding LBA and identify the use of inadequate outgroups as the reason previous deep genomic-scale analyses found strong support for Coelomata.},
  Doi                      = {10.1093/gbe/evq016},
  File                     = {:home/leo/articles/02-ImportedToMendeley/GenomeBiolEvol.2010.310.pdf:pdf},
  ISSN                     = {1759-6653},
  Keywords                 = {Animals,Bayes Theorem,Databases, Genetic,Eukaryota,Eukaryota: classification,Eukaryota: genetics,Evolution, Molecular,Genomics,Invertebrates,Invertebrates: classification,Invertebrates: genetics,Models, Genetic,Multigene Family,Phylogeny,Proteins,Proteins: classification,Proteins: genetics},
  Pmid                     = {20624736},
  Url                      = {http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2997542\&tool=pmcentrez\&rendertype=abstract}
}

# supermatrix, single gene approaches below:: 

@Article{Williams2012,
  Title                    = {{A congruent phylogenomic signal places eukaryotes within the Archaea.}},
  Author                   = {Williams, Tom A and Foster, Peter G and Nye, Tom M W and Cox, Cymon J and Embley, T Martin},
  Journal                  = {Proceedings of the The Royal Society / Biological Sciences},
  Year                     = {2012},
  Volume = {279},

  Month                    = oct,
  Number                   = {1749},
  Pages                    = {4870--4879},

  Abstract                 = {Determining the relationships among the major groups of cellular life is important for understanding the evolution of biological diversity, but is difficult given the enormous time spans involved. In the textbook 'three domains' tree based on informational genes, eukaryotes and Archaea share a common ancestor to the exclusion of Bacteria. However, some phylogenetic analyses of the same data have placed eukaryotes within the Archaea, as the nearest relatives of different archaeal lineages. We compared the support for these competing hypotheses using sophisticated phylogenetic methods and an improved sampling of archaeal biodiversity. We also employed both new and existing tests of phylogenetic congruence to explore the level of uncertainty and conflict in the data. Our analyses suggested that much of the observed incongruence is weakly supported or associated with poorly fitting evolutionary models. All of our phylogenetic analyses, whether on small subunit and large subunit ribosomal RNA or concatenated protein-coding genes, recovered a monophyletic group containing eukaryotes and the TACK archaeal superphylum comprising the Thaumarchaeota, Aigarchaeota, Crenarchaeota and Korarchaeota. Hence, while our results provide no support for the iconic three-domain tree of life, they are consistent with an extended eocyte hypothesis whereby vital components of the eukaryotic nuclear lineage originated from within the archaeal radiation.},
  Doi                      = {10.1098/rspb.2012.1795},
  File                     = {:home/leo/articles/02-ImportedToMendeley/ProcRSocB.2012.4870.pdf:pdf},
  ISSN                     = {1471-2954},
  Keywords                 = {00.UCA,algorithm,eukaryotes,evolution,gene trees,model selection,phylogenetics,tree of life,tree reconciliation},
  Mendeley-tags            = {00.UCA,algorithm,gene trees,model selection,tree reconciliation},
  Pmid                     = {23097517},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/23097517}
}
@Article{Lang2013,
  Title                    = {Phylogeny of Bacterial and Archaeal Genomes Using Conserved Genes: Supertrees and Supermatrices},
  Author                   = {Lang, Jenna Morgan and Darling, Aaron E. and Eisen, Jonathan A.},
  Journal                  = {PLoS ONE},
  Year                     = {2013},

  Month                    = {Apr},
  Number                   = {4},
  Pages                    = {e62510},
  Volume                   = {8},

  Doi                      = {10.1371/journal.pone.0062510},
  Editor                   = {Planet, Paul J.Editor},
  File                     = {:home/leo/articles/021/PLoSone.2013.0062510.pdf:PDF},
  Owner                    = {leo},
  Publisher                = {Public Library of Science},
  Timestamp                = {20131106.1809},
  Url                      = {http://dx.doi.org/10.1371/journal.pone.0062510}
}
@Article{Bininda2011,
  Title                    = {{Inferring the Tree of Life: chopping a phylogenomic problem down to size?}},
  Author                   = {Bininda-Emonds, Olaf R P},
  Journal                  = {BMC Biology},
  Year                     = {2011},

  Month                    = jan,
  Number                   = {Figure 2},
  Pages                    = {59},
  Volume                   = {9},

  Abstract                 = {The combination of molecular sequence data and bioinformatics has revolutionized phylogenetic inference over the past decade, vastly increasing the scope of the evolutionary trees that we are able to infer. A recent paper in BMC Biology describing a new phylogenomic pipeline to help automate the inference of evolutionary trees from public sequence databases provides another important tool in our efforts to derive the Tree of Life.},
  Doi                      = {10.1186/1741-7007-9-59},
  File                     = {:home/leo/articles/02-ImportedToMendeley/BMCBiol.2011.9-59.pdf:pdf},
  ISSN                     = {1741-7007},
  Owner                    = {leo},
  Pmid                     = {21933452},
  Timestamp                = {20130501.0000},
  Url                      = {http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3176480\&tool=pmcentrez\&rendertype=abstract}
}

@Article{Medina2011,
  Title                    = {{Reconstructing the fungal tree of life using phylogenomics and a preliminary investigation of the distribution of yeast prion-like proteins in the fungal kingdom.}},
  Author                   = {Medina, Edgar M and Jones, Gary W and Fitzpatrick, David A},
  Journal                  = {Journal of molecular evolution},
  Year                     = {2011},

  Month                    = oct,
  Number                   = {3-4},
  Pages                    = {116--33},
  Volume                   = {73},

  Abstract                 = {We have used three independent phylogenomic approaches (concatenated alignments, single-, and multi-gene supertrees) to reconstruct the fungal tree of life (FTOL) using publicly available fungal genomes. This is the first time multi-gene families have been used in fungal supertree reconstruction and permits us to use up to 66\% of the 1,001,217 genes in our fungal database. Our analyses show that different phylogenomic datasets derived from varying clustering criteria and alignment orientation do not have a major effect on phylogenomic supertree reconstruction. Overall the resultant phylogenomic trees are relatively congruent with one another and successfully recover the major fungal phyla, subphyla and classes. We find that where incongruences do occur, the inferences are usually poorly supported. Within the Ascomycota phylum, our phylogenies reconstruct monophyletic Saccharomycotina and Pezizomycotina subphyla clades and infer a sister group relationship between these to the exclusion of the Taphrinomycotina. Within the Pezizomycotina subphylum, all three phylogenies infer a sister group relationship between the Leotiomycetes and Sordariomycetes classes. However, there is conflict regarding the relationships with the Dothideomycetes and Eurotiomycetes classes. Within the Basidiomycota phylum, supertrees derived from single- and multi-gene families infer a sister group relationship between the Pucciniomycotina and Agaricomycotina subphyla while the concatenated phylogeny infers a poorly supported relationship between the Agaricomycotina and Ustilagomycotina. The reconstruction of a robust FTOL is important for future fungal comparative analyses. We illustrate this point by performing a preliminary investigation into the phyletic distribution of yeast prion-like proteins in the fungal kingdom.},
  Doi                      = {10.1007/s00239-011-9461-4},
  File                     = {:home/leo/articles/02-ImportedToMendeley/JMolEvol.2011.doi9461.pdf:pdf;:home/leo/articles/02-ImportedToMendeley/JMolEvol.2011.doi9461-4\_suppl1.pdf:pdf;:home/leo/articles/02-ImportedToMendeley/JMolEvol.2011.doi9461-4.pdf:pdf},
  ISSN                     = {1432-1432},
  Keywords                 = {Amino Acid Sequence,Bayes Theorem,Cluster Analysis,Computer Simulation,Fungal Proteins,Fungal Proteins: genetics,Fungi,Fungi: classification,Fungi: genetics,Genetic,Markov Chains,Models,Phylogeny,Prions,Prions: genetics,Sequence Alignment,fungal tree of life,markovian clustering \'{a},supertree \'{a} supermatrix \'{a},\'{a} fungal phylogeny \'{a}},
  Pmid                     = {21938499},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/21938499}
}
@Article{Romiguier2013,
  Title                    = {Less Is More in Mammalian Phylogenomics: AT-Rich Genes Minimize Tree Conflicts and Unravel the Root of Placental Mammals.},
  Author                   = {Romiguier, Jonathan and Ranwez, Vincent and Delsuc, Fr{\'{e}}d{\'{e}}ric and Galtier, Nicolas and Douzery, Emmanuel J P},
  Journal                  = {Mol. Biol. Evol.},
  Year                     = {2013},
  Number                   = {9},
  Pages                    = {2134-2144},
  Volume                   = {30},

  Doi                      = {10.1093/molbev/mst116},
  File                     = {file:///home/leo/articles/02-ImportedToMendeley/MolBiolEvol.2013.mst116.pdf:pdf},
  Keywords                 = {phylogenomics},
  Owner                    = {leo},
  Pmid                     = {23813978},
  Publisher                = {SMBE},
  Timestamp                = {2013.08.18},
  Url                      = {http://mbe.oxfordjournals.org/content/30/9/2134}
}

@Article{Salichos13,
  Title                    = {Inferring ancient divergences requires genes with strong phylogenetic signals.},
  Author                   = {Salichos, Leonidas and Rokas, Antonis},
  Journal                  = {Nature},
  Year                     = {2013},
  Number                   = {7449},
  Pages                    = {327-31},
  Volume                   = {497},

  Doi                      = {10.1038/nature12130},
  File                     = {:home/leo/articles/02-ImportedToMendeley/Nature.2013.nature12130.pdf:pdf;:home/leo/articles/02-ImportedToMendeley/Nature.2013.nature12130_suppl.pdf:pdf},
  Owner                    = {leo},
  Pmid                     = {23657258},
  Timestamp                = {2013.06.03_12:51}
}
@Article{Dunn2008,
  Title                    = {{Broad phylogenomic sampling improves resolution of the animal tree of life}},
  Author                   = {Dunn, Casey W and Hejnol, Andreas and Matus, D.Q. and Pang, K. and Browne, W.E. and Smith, S.A. and Seaver, E. and Rouse, G.W. and Obst, Matthias and Edgecombe, G.D. and Others},
  Journal                  = {Nature},
  Year                     = {2008},
  Number                   = {7188},
  Pages                    = {745--749},
  Volume                   = {452},

  File                     = {:home/leo/articles/02-ImportedToMendeley/Nature.2008.745.pdf:pdf},
  Publisher                = {Nature Publishing Group},
  Url                      = {http://www.nature.com/nature/journal/vaop/ncurrent/abs/nature06614.html}
}
@Article{Koonin2008,
  Title                    = {{Genomics of bacteria and archaea: the emerging dynamic view of the prokaryotic world.}},
  Author                   = {Koonin, Eugene V and Wolf, Yuri I},
  Journal                  = {Nucleic acids research},
  Year                     = {2008},

  Month                    = dec,
  Number                   = {21},
  Pages                    = {6688--719},
  Volume                   = {36},

  Abstract                 = {The first bacterial genome was sequenced in 1995, and the first archaeal genome in 1996. Soon after these breakthroughs, an exponential rate of genome sequencing was established, with a doubling time of approximately 20 months for bacteria and approximately 34 months for archaea. Comparative analysis of the hundreds of sequenced bacterial and dozens of archaeal genomes leads to several generalizations on the principles of genome organization and evolution. A crucial finding that enables functional characterization of the sequenced genomes and evolutionary reconstruction is that the majority of archaeal and bacterial genes have conserved orthologs in other, often, distant organisms. However, comparative genomics also shows that horizontal gene transfer (HGT) is a dominant force of prokaryotic evolution, along with the loss of genetic material resulting in genome contraction. A crucial component of the prokaryotic world is the mobilome, the enormous collection of viruses, plasmids and other selfish elements, which are in constant exchange with more stable chromosomes and serve as HGT vehicles. Thus, the prokaryotic genome space is a tightly connected, although compartmentalized, network, a novel notion that undermines the 'Tree of Life' model of evolution and requires a new conceptual framework and tools for the study of prokaryotic evolution.},
  Doi                      = {10.1093/nar/gkn668},
  File                     = {:home/leo/articles/02-ImportedToMendeley/NuclAcidsRes.2008.6688.pdf:pdf},
  ISSN                     = {1362-4962},
  Keywords                 = {Archaea,Archaea: metabolism,Bacteria,Bacteria: metabolism,Evolution, Molecular,Gene Transfer, Horizontal,Genes, Archaeal,Genes, Bacterial,Genome, Archaeal,Genome, Bacterial,Genomics,Signal Transduction},
  Pmid                     = {18948295},
  Url                      = {http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2588523\&tool=pmcentrez\&rendertype=abstract}
}
@Article{Sanderson2011,
  Title                    = {{Terraces in phylogenetic tree space.}},
  Author                   = {Sanderson, Michael J and McMahon, Michelle M and Steel, Mike},
  Journal                  = {Science (New York, N.Y.)},
  Year                     = {2011},

  Month                    = jul,
  Number                   = {6041},
  Pages                    = {448--50},
  Volume                   = {333},

  Abstract                 = {A key step in assembling the tree of life is the construction of species-rich phylogenies from multilocus--but often incomplete--sequence data sets. We describe previously unknown structure in the landscape of solutions to the tree reconstruction problem, comprising sometimes vast "terraces" of trees with identical quality, arranged on islands of phylogenetically similar trees. Phylogenetic ambiguity within a terrace can be characterized efficiently and then ameliorated by new algorithms for obtaining a terrace's maximum-agreement subtree or by identifying the smallest set of new targets for additional sequencing. Algorithms to find optimal trees or estimate Bayesian posterior tree distributions may need to navigate strategically in the neighborhood of large terraces in tree space.},
  Doi                      = {10.1126/science.1206357},
  File                     = {:home/leo/articles/02-ImportedToMendeley/Science.2011.411\_448.pdf:pdf},
  ISSN                     = {1095-9203},
  Keywords                 = {Algorithms,Angiosperms,Angiosperms: classification,Angiosperms: genetics,Animals,Arthropods,Arthropods: classification,Arthropods: genetics,Bayes Theorem,Biological Evolution,Colubridae,Colubridae: classification,Colubridae: genetics,Data Interpretation,Genomics,Genomics: methods,Likelihood Functions,Models,Nee,Phylogeny,Poaceae,Poaceae: classification,Poaceae: genetics,Sequence Alignment,Statistical},
  Mendeley-tags            = {Nee},
  Pmid                     = {21680810},
  Url                      = {http://www.ncbi.nlm.nih.gov/pubmed/21680810}
}

@Article{Huggins2011,
  Title                    = {{Bayes estimators for phylogenetic reconstruction.}},
  Author                   = {Huggins, P M and Li, W and Haws, D and Friedrich, T and Liu, J and Yoshida, R},
  Journal                  = {Systematic biology},
  Year                     = {2011},

  Month                    = jul,
  Number                   = {4},
  Pages                    = {528--40},
  Volume                   = {60},

  Abstract                 = {Tree reconstruction methods are often judged by their accuracy, measured by how close they get to the true tree. Yet, most reconstruction methods like maximum likelihood (ML) do not explicitly maximize this accuracy. To address this problem, we propose a Bayesian solution. Given tree samples, we propose finding the tree estimate that is closest on average to the samples. This "median" tree is known as the Bayes estimator (BE). The BE literally maximizes posterior expected accuracy, measured in terms of closeness (distance) to the true tree. We discuss a unified framework of BE trees, focusing especially on tree distances that are expressible as squared euclidean distances. Notable examples include Robinson-Foulds (RF) distance, quartet distance, and squared path difference. Using both simulated and real data, we show that BEs can be estimated in practice by hill-climbing. In our simulation, we find that BEs tend to be closer to the true tree, compared with ML and neighbor joining. In particular, the BE under squared path difference tends to perform well in terms of both path difference and RF distances.},
  Archiveprefix            = {arXiv},
  Arxivid                  = {0911.0645},
  Doi                      = {10.1093/sysbio/syr021},
  Eprint                   = {0911.0645},
  File                     = {:home/leo/articles/02-ImportedToMendeley/SystBiol.2011.528.pdf:pdf},
  ISSN                     = {1076-836X},
  Keywords                 = {Animals,Bayes Theorem,Computer Simulation,Genetic,Models,Phylogeny,Software,Urodela,Urodela: classification,Urodela: genetics},
  Owner                    = {leo},
  Pmid                     = {21471560},
  Timestamp                = {20131019.2115},
  Url                      = {http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3114872\&tool=pmcentrez\&rendertype=abstract http://www.ncbi.nlm.nih.gov/pubmed/21471560}
}
````
