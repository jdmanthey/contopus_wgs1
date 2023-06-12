# build scripts for stairway
java -cp stairway_plot_es Stairbuilder sordidulus.blueprint

java -cp stairway_plot_es Stairbuilder virens.blueprint

chmod +x sordidulus.blueprint.sh
./sordidulus.blueprint.sh

chmod +x virens.blueprint.sh
./virens.blueprint.sh
