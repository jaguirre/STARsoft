from setuptools import setup

__version__ = '0.0.1'

setup_args = {
    'name': 'kids',
    'author': 'James Aguirre',
    'author_email': 'jaguirre at sas.upenn.edu',
    'license': 'GPL',
    'package_dir' : {'kids' : 'src'},
    'packages' : ['kids'],
#    'scripts': glob.glob('scripts/*'),
    'version': __version__,
#    'package_data' : {'capo':['data/*.txt','data/psa32_apj/*.npz',\
#           'data/psa64_apj/*.npz','data/mwa128/*.dat',\
#           'data/21cmfast/ps*']}
}

if __name__ == '__main__':
#    from distutils.core import setup
    apply(setup, (), setup_args)
#    try:
#        from setuptools import setup
#    except:
#        from distutils.core import setup
#    try:
#        apply(setup, (), setup_args)
#    except:
#        setup(**setup_args)
