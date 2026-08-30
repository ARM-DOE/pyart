.. _overview:

Why Py-ART?
===========

The **Python ARM Radar Toolkit (Py-ART)** is a powerful and flexible open-source library designed for working with weather radar data. Whether you're a researcher, operational meteorologist, or data scientist, Py-ART provides the tools you need to analyze, visualize, and process radar data efficiently. Here's why Py-ART stands out:

---

A Common Data Model: The Radar Object
-------------------------------------

At the heart of Py-ART is its **common data model**, which is implemented through the `Radar` object. This unified data structure enables users to work with radar data from a variety of sources in a consistent and intuitive way. The benefits of this approach include:

- **Standardization Across Formats**: Radar data comes in many formats (e.g., NEXRAD, ODIM, UF, etc.), each with its own quirks and metadata conventions. Py-ART's `Radar` object abstracts these differences, providing a standardized interface for all supported formats.

- **Ease of Use**: By encapsulating radar data and metadata into a single object, Py-ART simplifies data manipulation and analysis. Users can access radar fields, coordinates, and metadata with minimal effort.

- **Interoperability**: The common data model facilitates integration with other scientific Python libraries (e.g., NumPy, SciPy, Matplotlib, xarray), enabling seamless workflows for data analysis and visualization.

- **Extensibility**: The `Radar` object is designed to be flexible, allowing users to add custom fields, attributes, or processing steps as needed.

---

xradar vs. Legacy Radar Data Structures
---------------------------------------

Py-ART represents a significant advancement over legacy radar data structures, particularly with its adoption of **xarray** through the `xradar` extension. Here's how Py-ART (with xradar) compares to older approaches:

.. list-table::
   :header-rows: 1
   :widths: 20 35 45

   * - **Feature**
     - **Legacy** (``pyart.io.read`` -> :class:`~pyart.core.Radar`)
     - **xradar-first** (``xd.io.open_*_datatree`` ->
       ``tree.pyart.to_radar()`` -> :class:`~pyart.xradar.Xradar`)
   * - Data Representation
     - Custom object with dictionaries
     - xarray/xradar ``DataTree`` under the hood, wrapped to duck-type
       ``Radar``
   * - Metadata Handling
     - Based on cfradial1 conventions
     - Based on cfradial2/FM301 conventions (via xradar)
   * - Format coverage
     - Py-ART's own readers (CfRadial, NEXRAD, Sigmet, MDV, ...)
     - Growing set of xradar readers; new format support is expected to
       land in xradar first, see :ref:`xradar_integration`
   * - Missing-data handling
     - ``numpy.ma.MaskedArray`` fields
     - Also ``numpy.ma.MaskedArray`` fields (NaN-masked on load)
   * - Interoperability
     - Minimal, package-by-package
     - Full integration with the PyData/xarray ecosystem

Not every :class:`~pyart.core.Radar` attribute/method has an
:class:`~pyart.xradar.Xradar` equivalent yet -- see the supported-API table
in :ref:`xradar_integration` for exactly what has been ported and verified.

The adoption of xarray in Py-ART allows for more intuitive and efficient handling of radar data, especially for large-scale or multi-dimensional datasets. It also aligns Py-ART with modern data science practices, ensuring long-term sustainability and usability.

---

Goals and Aspirations
---------------------

Py-ART was developed with the following goals and aspirations in mind:

1. **Accessibility**: Make radar data analysis accessible to a broad audience, from novice users to experienced researchers.
2. **Flexibility**: Provide tools that can handle diverse radar datasets and use cases, from basic visualization to advanced signal processing.
3. **Community-Driven Development**: Foster an open-source community where users can contribute new features, report issues, and collaborate on improvements.
4. **Scalability**: Enable efficient processing of large radar datasets, whether on a local machine or in a high-performance computing environment.
5. **Future-Proofing**: Stay aligned with modern data science tools and practices, ensuring Py-ART remains relevant as technology evolves.

By adhering to these principles, Py-ART aspires to be the go-to toolkit for radar data analysis, empowering users to unlock the full potential of their data.

---

In summary, Py-ART's common data model, modern architecture, and community-driven approach make it an indispensable tool for anyone working with radar data. Whether you're exploring new scientific questions or building operational workflows, Py-ART provides the foundation you need to succeed.
