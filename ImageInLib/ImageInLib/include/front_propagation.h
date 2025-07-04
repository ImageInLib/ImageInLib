
#pragma once

#include<iostream>
#include<vector>
#include "common_functions.h"

#define min_data(a,b)((a) < (b)) ? (a) : (b)
#define max_data(a,b)((a) > (b)) ? (a) : (b)

	typedef struct {
		size_t x, y;
		dataType arrival;
	}pointFastMarching2D;

	typedef struct {
		dataType K; //edge detection coef
		dataType thres;//edge detector threshold
		dataType eps; //path smothing parameter
		double radius;
	} Potential_Parameters;

	/// <summary>
	/// Swaps the values of two 2D points represented by pointFastMarching2D structures.
	/// </summary>
	/// <param name="a">Pointer to the first 2D point to swap.</param>
	/// <param name="b">Pointer to the second 2D point to swap.</param>
	void swap2dPoints(pointFastMarching2D* a, pointFastMarching2D* b);

	/// <summary>
	/// Restores the heap property by moving the element at the specified position down in a 2D fast marching heap.
	/// </summary>
	/// <param name="in_Process">A reference to the vector representing the heap of 2D fast marching points.</param>
	/// <param name="pos">The index of the element to heapify down.</param>
	void heapifyDown2D(std::vector<pointFastMarching2D>& in_Process, int pos);

	/// <summary>
	/// Restores the heap property by moving the element at the specified position up in a 2D fast marching heap.
	/// </summary>
	/// <param name="in_Process">A reference to the vector representing the heap of pointFastMarching2D elements.</param>
	/// <param name="pos">The index of the element to move up in the heap.</param>
	void heapifyUp2D(std::vector<pointFastMarching2D>& in_Process, int pos);

	/// <summary>
	/// Transforms a 2D vector of pointFastMarching2D objects into a heap structure.
	/// </summary>
	/// <param name="in_Process">A reference to a vector of pointFastMarching2D objects to be heapified.</param>
	void heapifyVector2D(std::vector<pointFastMarching2D>& in_Process);

	/// <summary>
	/// Deletes or processes the root element from a 2D fast marching heap represented by a vector of pointFastMarching2D objects.
	/// </summary>
	/// <param name="in_Process">A reference to a vector containing pointFastMarching2D objects that represent the heap.</param>
	void deleteRootHeap2D(std::vector<pointFastMarching2D>& in_Process);

	/// <summary>
	/// Adds a 2D point to the processing heap for fast marching algorithms.
	/// </summary>
	/// <param name="in_Process">A reference to the vector representing the current processing heap of 2D points.</param>
	/// <param name="point">The 2D point to be added to the processing heap.</param>
	void addPointHeap2D(std::vector<pointFastMarching2D>& in_Process, pointFastMarching2D point);
	
	/// <summary>
	/// Returns the index of the element at the specified 2D coordinates in a 1D vector representing a 2D heap structure.
	/// </summary>
	/// <param name="in_Process">A reference to a vector of pointFastMarching2D objects representing the 2D heap.</param>
	/// <param name="i">The row index of the element.</param>
	/// <param name="j">The column index of the element.</param>
	/// <returns>The index in the vector corresponding to the element at position (i, j) in the 2D heap.</returns>
	int getIndexFromHeap2D(std::vector<pointFastMarching2D>& in_Process, size_t i, size_t j);

	/// <summary>
	/// Return the minimum action value in x direction for a given point in a 2D action map.
	/// </summary>
	/// <param name="actionMapPtr">Pointer to the action map</param>
	/// <param name="length">image length</param>
	/// <param name="width">image width</param>
	/// <param name="ind_x">x coordinate of current point</param>
	/// <param name="ind_y">y coordinate of current point</param>
	/// <returns>INIFITY if the point is on the boundary or UNKNOW or minimum value using upwind principle</returns>
	dataType selectX(dataType* actionMapPtr, const size_t length, const size_t width, const size_t ind_x, const size_t ind_y);

	/// <summary>
	/// Return the minimum action value in y direction for a given point in a 2D action map.
	/// </summary>
	/// <param name="actionMapPtr">Pointer to the action map</param>
	/// <param name="length">image length</param>
	/// <param name="width">image width</param>
	/// <param name="ind_x">x coordinate of current point</param>
	/// <param name="ind_y">y coordinate of current point</param>
	/// <returns>INIFITY if the point is on the boundary or UNKNOW or minimum value using upwind principle</returns>
	dataType selectY(dataType* actionMapPtr, const size_t length, const size_t width, const size_t ind_x, const size_t ind_y);

	dataType solve2dQuadratic(dataType X, dataType Y, dataType P, PixelSpacing h);

	bool frontPropagation2D(Image_Data2D imageData, dataType* distancePtr, dataType* potentialPtr, Point2D* seedPoints);
