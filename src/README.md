# linearham                          {#mainpage}

## Handy links

* [linearham GitHub repository](https://github.com/matsengrp/linearham)
* [Catch tutorial](https://github.com/philsquared/Catch)
* [Eigen docs](http://eigen.tuxfamily.org/dox/group__QuickRefPage.html)
* [Google C++ style guide](https://google.github.io/styleguide/cppguide.html)


## const, references, and Refs
### Function arguments
All passing of Eigen objects into functions should be by reference to avoid copying.
There is an [Eigen construct called `Ref`](http://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html#TopicUsingRefClass) that allows passing either a matrix-like object or a `Block` of that object.
There are a silly number of combinations of these things, but we will use the following two for input (see [this SO question](http://stackoverflow.com/questions/21132538/correct-usage-of-the-eigenref-class), and the examples in [the Ref docs](https://eigen.tuxfamily.org/dox-devel/classEigen_1_1Ref.html)):

* For const input: `const Eigen::Ref<const Eigen::MatrixXd>&`
* For read/write input: `Eigen::Ref<Eigen::MatrixXd>`

### Returns
* If an Eigen object is made and then returned, it won't be copied because of [return value optimization](https://en.wikipedia.org/wiki/Return_value_optimization).
* For a const accessor: `const Eigen::MatrixXd&`
* For now we don't need accessors that can be used as an lvalue (otherwise why not just make the member variable public?), but I think we would return an `Eigen::Ref<Eigen::MatrixXd>`.

Note that trying to return a  `const Eigen::Ref<const Eigen::MatrixXd>&` as an accessor, the compiler gives an error saying that we can't return a reference to a local temporary object, which I think means that we'd have to make a `Ref` object on the fly when what we really have is a Matrix or something.



## Glossary

* *segment:* a contiguous sequence of bases that come from one type of state, e.g. V, N, D, etc.
* *symbol:* what gets emitted. Here nucleotide bases.


## Smooshables

The fundamental abstraction in linearham is a *smooshable*.
A smooshable is a linear segment of sequence that has a probability associated with each pair of begin and end points.
For example, a smooshable coming from a read aligned to a given germline sequence would have probability for a given pair of begin and end points equal to the probability of emitting the corresponding segment of read from that germline gene.
That is, it would be the probability of entering the linear segment of the HMM at the specified place, emitting the observed bases, and exiting at the specified place.

The allowable range for the begin and end points are called the *left flex* and *right flex*, respectively.
The begin and end points are indexed from the beginning of the smooshable.
An end point of zero corresponds to `right_flex` from the beginning, like so:

![](http://i.imgur.com/7CKeqed.png)

To *smoosh* two smooshables is to juxtapose two of them and consider all the options of going from the left one to the right one.
Here "consider" means both calculating marginal probabilities (summing over transition points) and calculating Viterbi paths.
