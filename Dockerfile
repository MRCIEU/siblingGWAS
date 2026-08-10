FROM rocker/r-ver:4.5.2

RUN Rscript -e 'install.packages( \
        c( \
            "data.table", \
            "lmtest", \
            "sandwich", \
            "R.utils" \
        ), \
        repos="https://www.stats.bris.ac.uk/R/")'

COPY within-family.R /scripts/
COPY within-family-plots.R /scripts/

COPY scripts/within-family.sh /usr/local/bin/within-family
COPY scripts/within-family-plots.sh /usr/local/bin/within-family-plots

RUN chmod +r /scripts/within-family.R /scripts/within-family-plots.R && \
    chmod +x /usr/local/bin/within-family /usr/local/bin/within-family-plots

CMD ["within-family", "-h"]
