#!/bin/bash
cd data

if [ -d "db-light" ]; then
    echo "Database exists. Skip download."
else
    echo "Started loading bakta database"
    wget https://zenodo.org/record/7669534/files/db-light.tar.gz
    tar -xzf db-light.tar.gz
    rm db-light.tar.gz
fi

./amrfinder_update --force_update --database db-light/amrfinderplus-db/
