#!/bin/bash

export AZ_BLOB_PREFIX=snakemake # This must be the container name
export AZ_BLOB_ACCOUNT_URL="https://${stgacct}.blob.core.windows.net/$sastoken" # Get this information from the Azure Portal or Azure CLI

snakemake --config use_azure=true --kubernetes \
    --default-remote-prefix $AZ_BLOB_PREFIX --default-remote-provider AzBlob \
    --envvars AZ_BLOB_ACCOUNT_URL AZ_BLOB_PREFIX --use-conda --jobs 2 --verbose