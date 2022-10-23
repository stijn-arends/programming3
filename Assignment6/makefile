# Makefile for the project

# Links:
# manage python envs using makefiles: https://ljvmiranda921.github.io/notebook/2021/05/12/how-to-manage-python-envs/
# Makefile tutorial for python: https://medium.com/analytics-vidhya/makefile-make-rules-and-automate-in-python-982a41c99275
# Activate venv inside makefile: https://stackoverflow.com/questions/33839018/activate-virtualenv-in-makefile 
# More advanced example: https://dev.to/mmphego/why-you-should-add-makefile-into-your-python-project-20j2

.SILENT: dev requirements create_venv

# Does not work at the moment, idea was to silence the stdout of the commands
# stdout=false

# conditional:
# ifeq ($(stdout), false)
# 	echo "Hello there"
# 	@redirect=2> /dev/null 1>&2
# else
# 	echo "Ow no"
# 	redirect= 
# endif

# If you want to silence the stdout and stderr of commands add this to the end of the command: 2> /dev/null 1>&2

# Variables
virtual_env=~/.venv/pubmed/
install=install/

create_venv:
	@python -m venv $(virtual_env) 2> /dev/null 1>&2
	@$(virtual_env)bin/python -m pip install --upgrade pip 2> /dev/null 1>&2
	@$(virtual_env)bin/pip3 install pip-tools 2> /dev/null 1>&2

requirements: create_venv $(install)requirements.in ## create dependencies for release
	@$(virtual_env)bin/pip-compile -o $(install)requirements.txt $(install)requirements.in 2> /dev/null 1>&2

dev: requirements create_venv # creates a development environment, install deps
	@$(virtual_env)bin/pip-sync $(install)requirements.txt 2> /dev/null 1>&2


venv: dev # to call python script add: && python script.py
	. $(virtual_env)/bin/activate 