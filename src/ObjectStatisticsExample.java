import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import annotool.ImgDimension;
import annotool.extract.FeatureExtractor;
import annotool.io.DataInput;

/*
 * Statistics regarding objects in the image
 * 
 * 1. The number of objects
 * 2. The average number of above-threshold pixel per object
 * 3. The variance of the number of above-threshold pixels per object
 * 4. The ratio of the size of the largest object to the smallest.
 * 5. The average object distance to the CenterOfMass (CoM)
 * 6. The variance of the object distance from the image CoM.
 * 7. The ratio of the largest to the smallest object to image CoM.
 * 
 * Ref: "A neural network classifier capable of recognizing the patterns of all major subcellular
 *    structures in fluorescence microscope images of HeLa cells",  M.V. Boland and R. F. Murphy
 *    Bioinformatics, 2001.
 */
public class ObjectStatisticsExample implements FeatureExtractor {

	//properties about the incoming data
	int length;
	int imageType;   
	DataInput problem = null;
	ArrayList allData = null;
	//size of the current image being processed
	int width, height; 

	//features to return
	protected float[][] features = null;

	//masks for sign extension
	final int mask = 0xff; 
	final int longmask = 0xffff;
	final int NUMBER_OF_FEATURES = 7;

	public void setParameters(HashMap<String, String> para) {
		// Use if needed

	}

	public float[][] calcFeatures(DataInput problem) throws Exception {

		this.problem = problem;
		this.length = problem.getLength();
		if (problem.ofSameSize() != false)
		{
			this.width  =  problem.getWidth();
			this.height  = problem.getHeight();
		}
		this.imageType = problem.getImageType();
		features = new float[length][NUMBER_OF_FEATURES];

		return calcFeatures();
	}

	public float[][] calcFeatures(ArrayList data, int imageType,
			ImgDimension dim) throws Exception 
			{
		// To be done.
		throw new Exception("To be implemented later.");
			}

	public boolean is3DExtractor() {
		return false;
	}

	//calculate and fill the result array
	private float[][] calcFeatures() throws Exception
	{
		Object data; //The 2D image. It can be byte[] or short[] for float[] depending on image type
		int[] tag;   //tags of found objects. See getObjectsInImage() for detail.
		HashMap<Integer, Integer> tagCount;  //the size in # of pixels for each object 
		java.util.Set<Integer> keySet;

		for (int image_num = 0; image_num < this.length; image_num++) {
			if (image_num % 50 == 0) System.out.println();
			System.out.print(image_num+" ");

			//get dimension for current image
			if (problem !=null)
			{	
				if (!problem.ofSameSize())
				{  	//set the size for this image
					this.height = problem.getHeightList()[image_num];
					this.width = problem.getHeightList()[image_num];
				}
			}
			//get data -- current image
			if (problem !=null)
				data = problem.getData(image_num,1);
				//data = problem.getData(image_num);
			else  //alldata is not null
				data = allData.get(image_num);

			//get threshold for binarizing the image.  You can use your own way here if you want.
			ij.ImagePlus iplus = problem.getImagePlus(image_num);
			int threshold = iplus.getProcessor().getAutoThreshold();

			//get objects for current image
			tag = getObjectsInImage(data, width, height, threshold);
			tagCount = getTagCount(tag, width, height);
			keySet = tagCount.keySet();
			int[] center = getCofM(data, width, height);

			//set the features
			features[image_num][0] =  calcNumberOfObjects(tagCount);
			features[image_num][1] =  calcAvgSizeOfObjects(tagCount);
			features[image_num][2] =  calcVarSizeOfObjects(tagCount, features[image_num][1]);
			features[image_num][3] =  calcRatioBigToSmallOfObjects(tagCount);
			features[image_num][4] =  calcAvgDistanceToCoM(data, width, height, keySet, tag, center);
			features[image_num][5] =  calcVarDistanceToCoM(data, width, height, keySet, tag, center, features[image_num][4]);
			features[image_num][6] =  calcRatioFurthestToClosestToCoM(data, width, height, keySet, tag, center);
		}

		return features;
	}

	//**********************   METHODS  TO BE FILLED BY YOU   *****************************************/

	private float calcNumberOfObjects(HashMap<Integer, Integer>tagCount)
	{
		return tagCount.size();
	}


	private float calcAvgSizeOfObjects(HashMap<Integer, Integer> tagCount)
	{
		float avgSize=0;
		for(Entry<Integer, Integer> entry : tagCount.entrySet()){
			avgSize += entry.getValue();
		}
		avgSize /= tagCount.size();
		return avgSize;	
	}

	private float calcVarSizeOfObjects(HashMap<Integer, Integer> tagCount, float avg)
	{
		float variance=0;
		for(Entry<Integer, Integer> entry : tagCount.entrySet()){
			variance+=(avg-entry.getValue())*(avg-entry.getValue());
		}
		return variance/tagCount.size();
	}

	private float calcVarSizeOfObjects(HashMap<Integer, Integer> tagCount)
	{
		float avgSize=0;
		for(Entry<Integer, Integer> entry : tagCount.entrySet()){
			avgSize += entry.getValue();
		}
		avgSize /= tagCount.size();
		float variance=0;
		for(Entry<Integer, Integer> entry : tagCount.entrySet()){
			variance+=(avgSize-entry.getValue())*(avgSize-entry.getValue());
		}
		return variance/tagCount.size();
	}

	private float calcRatioBigToSmallOfObjects(HashMap<Integer, Integer> tagCount)
	{
		int max = Integer.MIN_VALUE;
		int min = Integer.MAX_VALUE;
		for(Entry<Integer, Integer> entry : tagCount.entrySet()){
			if(entry.getValue() < min)
				min = entry.getValue();
			if (entry.getValue() > max)
				max = entry.getValue();
		}
		return max/min;
	}

	private float calcAvgDistanceToCoM(Object imgData, int width, int height, java.util.Set<Integer> keySet, int[] tag, int[] center)
	{
		
		return 0;
	}

	private float calcVarDistanceToCoM(Object imgData, int width, int height, java.util.Set<Integer> keySet, int[] tag, int[] center, float avgDist)
	{
		return 0;
	}

	private float calcRatioFurthestToClosestToCoM(Object imgData, int width, int height, java.util.Set<Integer> keySet, int[] tag, int[] center)
	{
		return 0;		
	}

	//**********************   HELPER METHODS   *****************************************/

	//this method go through the images and get the objects of interest
	//It  stores each object's information for later processing, in an array with tag
	// threshold is the binarization threshold.
	//The basic algorithm is similar as ImageJ plugin 3D Object Counter
	//http://rsbweb.nih.gov/ij/plugins/track/objects.html
	private int[]  getObjectsInImage(Object imgData, int width, int height, int threshold)
	{
		int x, y;
		int value;
		int arrayIndex;
		int[] tag = new int[width*height];
		boolean[] isObjPixel = new boolean[width*height];

		//binarize
		arrayIndex = 0;
		for(y = 0; y < height; y++) {
			for(x = 0; x < width; x++) {
				//get value
				value = getValue(imgData, arrayIndex);
				if(value < threshold)
					isObjPixel[arrayIndex] = true;
				else
					isObjPixel[arrayIndex] = false;
				arrayIndex++;
			}
		}

		int tagvois;
		int ID = 1;
		int minTag;
		int i, offset;
		int nX = -1, nY = -1;
		arrayIndex = 0;

		//First ID attribution
		for(y = 0; y < height; y++) {
			for(x = 0; x < width; x++) {
				if(isObjPixel[arrayIndex]) {
					tag[arrayIndex] = ID;
					minTag = ID;

					i = 0;
					//Find the minimum tag in the neighbors pixels
					for (nY = y - 1; nY <= y + 1; nY++) {
						for (nX = x - 1; nX <= x + 1; nX++) {
							if(withinBounds(nX, nY, width, height)) {
								offset = offset(nX, nY, width);
								if(isObjPixel[offset]) {
									i++;	//If neighbor pixel is object, increment neighbor object pixel count
									tagvois = tag[offset];
									if (tagvois != 0 && tagvois < minTag) 	//If any neighbor object pixel is already tagged, use that
										minTag = tagvois;
								}
							}
						}
					}

					tag[arrayIndex] = minTag;
					if (minTag == ID)
					{  //If current pixel was given the new tag, increment tag number for next run
						ID++;
					}
				}
				arrayIndex++;
			}
		}

		//Minimization of IDs = connection of structures
		arrayIndex = 0;
		for(y = 0; y < height; y++) 
		{
			for(x = 0; x < width; x++) 
			{
				if(isObjPixel[arrayIndex]) 
				{
					minTag = tag[arrayIndex];

					//Find the minimum tag in the neighbors pixels
					for (nY = y - 1; nY <= y + 1; nY++) {
						for (nX = x - 1; nX <= x + 1; nX++) {
							if(withinBounds(nX, nY, width, height)) {
								offset = offset(nX, nY, width);
								if(isObjPixel[offset]) {
									tagvois = tag[offset];
									if (tagvois != 0 && tagvois < minTag)
										minTag = tagvois;
								}
							}
						}
					}

					//Replacing tag by the minimum tag found
					for (nY = y - 1; nY <= y + 1; nY++) {
						for (nX = x - 1; nX <= x + 1; nX++) {
							if(withinBounds(nX, nY, width, height)) {
								offset = offset(nX, nY, width);
								if(isObjPixel[offset]) {
									tagvois = tag[offset];
									if (tagvois != 0 && tagvois != minTag)
										replacetag(tag, tagvois,minTag);
								}
							}
						}
					}
				}
				arrayIndex++;
			}
		}

		return tag;
	}

	//get pixel value
	private int getValue(Object imgData, int arrayIndex)
	{
		int value =0;
		if (imageType == DataInput.GRAY8 || imageType ==  DataInput.COLOR_RGB)
			value = (int) (((byte[])imgData)[arrayIndex] & mask);
		else if (imageType == DataInput.GRAY16)
			value = (int) (((short[])imgData)[arrayIndex] & longmask);
		return value;
	}

	private boolean withinBounds(int x, int y, int width, int height) {
		return (x >= 0 && x < width && y >= 0 && y < height);
	}

	//Calculates the 1D array index for given x and y
	private int offset(int x, int y, int width) {
		return x + y * width;
	}

	public void replacetag(int[] tag, int m, int n){
		for (int i=0; i < tag.length; i++) 
			if (tag[i] == m) tag[i] = n;
	}

	//build a hashmap that stores the number of pixels for each object tag.	
	private HashMap<Integer, Integer> getTagCount(int[] tag, int width, int height) {
		//Key: object ID, value: pixel count
		HashMap<Integer, Integer> tagCount = new HashMap<Integer, Integer>();

		int arrayIndex = 0;
		for(int y=0; y < height; y++) {
			for(int x = 0; x < width; x++) {
				if(tag[arrayIndex] != 0) {
					int count = 0;
					if(tagCount.containsKey(tag[arrayIndex])) {
						count = tagCount.get(tag[arrayIndex]);
					}
					count++;
					tagCount.put(tag[arrayIndex], count);
				}

				arrayIndex++;
			}
		}

		return tagCount;
	}

	//get the coordinates of center of mass
	private int[] getCofM(Object imgData, int width, int height)
	{
		int arrayIndex = 0;
		int sumx = 0, sumvaluex = 0;
		int sumy = 0, sumvaluey = 0;
		int value;
		for(int y = 0; y < height; y++) {
			for(int x = 0; x < width; x++) {
				value = getValue(imgData, arrayIndex);
				sumvaluex +=value*x;
				sumvaluey +=value*y;
				sumx += value;
				sumy += value;
				arrayIndex++;
			}
		}
		int[] center = new int[2];
		center[0] = sumvaluex/sumx;
		center[1] = sumvaluey/sumy;
		return center;
	}

	//get CofM for a particular object
	private int[] getCofM(Object imgData, int width, int height, int[] tag, int index)
	{
		int arrayIndex = 0;
		int sumx = 0, sumvaluex = 0;
		int sumy = 0, sumvaluey = 0;
		int value;
		for(int y = 0; y < height; y++) 
		{
			for(int x = 0; x < width; x++) 
			{
				if(tag[arrayIndex] == index)
				{
					value = getValue(imgData, arrayIndex);
					sumvaluex +=value*x;
					sumvaluey +=value*y;
					sumx += value;
					sumy += value;
				}
				arrayIndex++;
			}
		}

		int[] center = new int[2];
		center[0] = sumvaluex/sumx;
		center[1] = sumvaluey/sumy;
		return center;
	}

	public static void main(String[] args)
	{
		//to be changed by you 
		String directory="/Users/PARTHASARATHY/Desktop/PatterRec/FinalProjectPattern/src/data/";
		String ext = ".png";
		String ch = "g";	

		try
		{
			annotool.io.DataInput problem = new annotool.io.DataInput(directory, ext, ch, true);
			FeatureExtractor fe = new ObjectStatisticsExample();
			//print out
			float[][] features = fe.calcFeatures(problem);
			System.out.println();
			for(int i = 0; i < features.length; i++)
			{
				System.out.print("Image #" + i + ": ");
				for(int j=0; j< features[i].length; j++)
					System.out.print(features[i][j]+" ");
				System.out.println();
			}
		}catch(Exception e)
		{
			e.printStackTrace();
			System.out.println(e.getMessage());
		}
	}

}