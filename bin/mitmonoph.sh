if [ -z $CMSSW_BASE ]
then
  echo ""
  echo " Setting up MitHgg failed! (\$CMSSW_BASE = is empty)."
  echo ""
else
  export MIT_MONOPH_DIR="$CMSSW_BASE/src/MitMonoPhoton"
  export PATH="$MIT_MONOPH_DIR/bin:${PATH}"
  export PYTHONPATH="$MIT_MONOPH_DIR/python:${PYTHONPATH}"
fi

