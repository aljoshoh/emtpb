#!/usr/bin/env bash

echo "Authenticating at: $(date)" | tee -a ~/rstudio_auth_r4.log
echo "$# arguments given: $@" | tee -a ~/rstudio_auth_r4.log

# Confirm username is supplied
if [ -z "$1" ]
  then
    echo "Argument 1 is empty, it should contain USERNAME" | tee -a ~/rstudio_auth_r4.log
fi
USERNAME="${1}"
echo "Authenticating user ${USERNAME}" | tee -a ~/rstudio_auth_r4.log

# Confirm user environment variable exists
if [[ -z ${RSTUDIO_USER} ]]; then
  echo "The environment variable RSTUDIO_USER is not set" | tee -a ~/rstudio_auth_r4.log
  exit 1
fi

# Confirm password environment variable exists
if [[ -z ${RSTUDIO_PASSWORD} ]]; then
  echo "The environment variable RSTUDIO_PASSWORD is not set" | tee -a ~/rstudio_auth_r4.log
  exit 1
fi

# Read in the password from user
read -s -p "Password: " PASSWORD
echo ""

if [[ ${USERNAME} == ${RSTUDIO_USER} && ${PASSWORD} == ${RSTUDIO_PASSWORD} ]]; then
  echo "Authentication successful!" | tee -a ~/rstudio_auth_r4.log
  exit 0
else
  echo "Invalid authentication" | tee -a ~/rstudio_auth_r4.log
  exit 1
fi
