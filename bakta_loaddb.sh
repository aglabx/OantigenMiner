#!/bin/bash
cd data

if [ -d "db-light" ]; then
    echo "Database exists. Skip download."
else
    echo "Started loading bakta database"
    wget https://zenodo.org/record/7669534/files/db-light.tar.gz
    tar -xzf db-light.tar.gz
    # rm db-light.tar.gz
fi

if [ -d "db-light/amrfinderplus-db" ]; then
    echo amrfinderplus-db exists too
else
    amrfinder_update --force_update --database db-light/amrfinderplus-db/
fi