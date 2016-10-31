# abcngspipelines

### Who do I talk to? ###

* Alexis Blanchet-Cohen
* alexis.blanchet-cohen@mail.mcgill.ca

The IRCM NGS pipelines are a collection of Python and R scripts for the analysis of next-generation sequencing data.
The scripts form a flexible and modular framework for next generation sequencing analyses.
The pipelines currently run on one of Compute Canada's computing clusters, Guillimin, and on IRCM's local server, Gen01, but they can be deployed on other servers. 

The scripts are still under active development.
As such, they are not quite ready for release to the general public.
They are, however, already used to handle the large amount of next-generation sequencing analyses carried out by the IRCM Bioinformatics Core Facility.
They are also used to train newly arrived bioinformaticians at IRCM. 

The Python scripts simply generate the command scripts, but do not submit commands directly to the queue.
This allows an extra level of control, since the command scripts can be edited directly.

### Full documentation (not up-to-date) ###
* http://confluence.ircm.qc.ca/display/IP/IRCM+NGS+pipelines

