from setuptools import setup

setup(
    name='vcaller',
    version='0.1',
    py_modules=['vcaller'],
    include_package_data=True,
    install_requires=[
        'click',
    ],
    entry_points='''
        [console_scripts]
        vcaller=vcaller:cli
    ''',
)