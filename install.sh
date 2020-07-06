if [ "$1" != "" ]; then
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	chmod +x Miniconda3-latest-Linux-x86_64.sh
	./Miniconda3-latest-Linux-x86_64.sh -b -p $1/miniconda3
	source $1/miniconda3/etc/profile.d/conda.sh
	conda env create -n procare -f procare_environment.yml
	conda activate procare
	pip install procare_python_package/

	echo "source $1/miniconda3/etc/profile.d/conda.sh" > activate.sh
	echo "conda activate procare" >> activate.sh
	#chmod +x activate.sh
else
	echo "please provide install directory. For example bash install.sh \$HOME"
fi