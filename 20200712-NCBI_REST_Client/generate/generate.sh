# RA, 2020-07-12

# Auto-generate a python client for the NCBI REST
# https://api.ncbi.nlm.nih.gov/datasets/v1alpha/

# Note:
# As of 2020-07-12 there is this auto-generated client
# https://github.com/ncbi/datasets/tree/master/client_docs/python

rm -r client

wget -nc https://repo1.maven.org/maven2/io/swagger/codegen/v3/swagger-codegen-cli/3.0.20/swagger-codegen-cli-3.0.20.jar -O swagger-codegen-cli.jar
wget -nc https://api.ncbi.nlm.nih.gov/datasets/v1alpha/datasets.openapi.yaml -O specs.yaml
java -jar swagger-codegen-cli.jar generate \
  -l python \
  -i specs.yaml \
  -o client \
  -c config.json
