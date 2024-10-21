# Erratum for Bowyer et al., 2019 DOI: 10.3390/microorganisms7010017

This git repo contains what we had hoped would be a published erratum for the article published in *Microorganisms*: ‘Socioeconomic Status and the Gut Microbiome: A TwinsUK Cohort’ (Microorganisms, 2019 Jan; 7(1): 17. DOI: 10.3390/microorganisms7010017). We have contacted the journal editors several times over the last few years since we noted the error, but unfortunately do not believe myself or the corresponding author have received a reply, and so are sharing publicly here. Please do contact me at ruth.c.bowyer [at] kcl.ac.uk if we have missed this opportunity to correct within the jounral, or there is another pathway to doing so! 

The error arose from the misspecification of variable order within a model; the corrective note is below, and an updated script can be found in ```/Updated Script``` above, comparing both the incorrect and correct output for the first model, and the correctly specified outputs for the remaining as an updated ‘additional file 4’. The script we hope is particularly of use as we know that other researchers make use of the scripts we published with the original article. 

**Importantly, the error does not affect the conclusions of the paper, and we do not directly discuss the erroneous result in the discussion or abstract sections**. Regardless, we sincerely apologise for the error and share here in the sincere belief of scientific accountability and honesty. 

## Corrective Note 

In the article ‘Socioeconomic Status and the Gut Microbiome: A TwinsUK Cohort Study’ by Ruth C E Bowyer, Matthew A Jackson, Caroline I Le Roy, Mary Ni Lochlainn, Tim D Spector, Jennifer B Dowd and Claire J Steves (Microorganisms, 2019 Jan; 7(1): 17. DOI: 10.3390/microorganisms7010017) there was a misspecification of a model relating to the adjusted results reported as page 5 as follows:

> Examining intra-individual (beta) microbiome diversity, we found significant differences across education and IMD groups in crude and adjusted NPMANOVA for Bray–Curtis dissimilarity (Education: crude F(3, 1425) = 1.71, p = 0.0004; adjusted F(3, 1425) = 1.75, p = 0.0004, IMD: crude F(4, 1671) = 1.37, p = 0.008; adjusted F(4, 1671) = 1.41, p = 0.006, permutations = 5000) and weighted UniFrac (Education: crude F(3, 1425) = 1.74, p = 0.0022; adjusted F(3, 1425) = 1.73, p = 0.0028, IMD: crude F(4, 1671) = 1.35, p = 0.03, adjusted F(4, 1671) = 1.38, p = 0.02).

The order of variables added to the model was incorrect, in that the primary variable of interest (education/IMD/income) should have been specified last. On re-running the models in the correct order, only adjusted models for education were significant and thus the paragraph should instead read: 

> Examining intra-individual (beta) microbiome diversity, we found significant differences between IMD groups in crude models, and in education groups for crude and adjusted NPMANOVA for Bray–Curtis dissimilarity (Education: crude F(3, 1425) = 1.71, p = 0.0004; adjusted F(3, 1425) = 1.48, p = 0.007, IMD: crude F(4, 1671) = 1.37, p = 0.008; permutations = 5000) and weighted UniFrac (Education: crude F(3, 1425) = 1.74, p = 0.0022; adjusted F(3, 1425) = 1.38, p = 0.043,  IMD: crude F(4, 1671) = 1.35, p = 0.03). 

An updated script for correct specification of the adjusted models is attached. The rest of the article remains unchanged. 
