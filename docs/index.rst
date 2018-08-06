.. datamaterials_neighbors documentation master file, created by
   sphinx-quickstart on Sat Aug  4 18:54:09 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to datamaterials_neighbors's documentation!
===================================================

Structure
---------

.. automodule:: datamaterials_neighbors.structure
   
   .. rubric:: Primary methods

   .. autosummary::

      from_file
      from_parameters

   .. autosummary::
   
   .. rubric:: Functions

   .. autosummary::
   
      get_subspecies_labels
      label_subspecies
   
   .. rubric:: Classes

   .. autosummary::
   
      StructureParameters
   

Neighbors
---------

.. automodule:: datamaterials_neighbors.neighbors
  
   .. rubric:: Primary method

   .. autosummary::
   
      count_neighbors

   .. rubric:: Functions

   .. autosummary::
   
      add_subspecie_labels_if_missing
      append_site_i_neighbor_distance_data
      count_neighbors_within_distance_groups
      define_bin_intervals
      define_bins_to_group_and_sort_by_distance
      extract_neighbor_distance_data
      find_unique_distances
      get_neighbor_distances_data_frame
      group_site_index_pairs_by_distance
   
   .. rubric:: Classes

   .. autosummary::
   
      NeighborData


Interactions
------------

.. automodule:: datamaterials_neighbors.interactions

   .. rubric:: Primary method

   .. autosummary::
   
      build_model

   .. rubric:: Functions

   .. autosummary::
   
      aggregate_interaction_coefficients
      build_magnetic_patterns_data_frame
      compute_interaction_signs
      compute_model_coefficients
      group_subspecie_pairs_and_rank_by_distance
      label_interaction_parameters
      multiply_interaction_signs_and_neighbor_count
      spread_parameter_name_column


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
